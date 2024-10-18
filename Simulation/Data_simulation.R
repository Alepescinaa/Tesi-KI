################################
# Upload libraries & functions
################################

knitr::opts_chunk$set(echo = TRUE)

library(dplyr)
library(haven)
library(mstate)
library(flexsurv)
library(lubridate)
library(hesim)
library(data.table)
library(ggplot2)
library(future)
library(future.apply)
library(parallel)

source("./functions/prepare_data.R")
source("./functions/simulation.R")
source("./functions/apply_right_censoring.R")
source("./functions/interval_censoring.R")
source("./functions/run_simulation.R")
source("./functions/interval_censoring_id.R")

################################################
# Upload data and preprocess them in long format
################################################

load("SNAC_data_Alessandra.RData")
data_raw <- SNAC_K_dataset

data_raw <- prepare_data(data_raw) 


# uptowave3 <- data_raw[complete.cases(data_raw[, c("Date_wave1", 
#                                                       "Date_wave2", 
#                                                       "Date_wave3")]), ]
# uptowave3 <- uptowave3[uptowave3$age_wave1>78,]
# 
# diff1 <- as.numeric(uptowave3$Date_wave2 - uptowave3$Date_wave1)/365
# diff2 <- as.numeric(uptowave3$Date_wave3 - uptowave3$Date_wave2)/365
# tot <- (diff1+diff2)/2
# mean(tot)
# 
# temp <- data.frame(uptowave3$lopnr, diff1, diff2)

tmat <- mstate::transMat(x = list(c(2, 3),c(3),c()), names = c("Dementia-free","Dementia", "Death")) 
print(tmat)

#converting smoking to bivariate covariate
# 0 has never smoked, 1 used to smoke, 2 current smoker
# we set 0 doesn't smoke, 1 current smoker
data_raw$Smoking <- if_else(data_raw$Smoking==0|data_raw$Smoking==1,0,1) 

#converting education to bivariate covariate
# 0 elementary, 1 high school, 2 uni or more
# we set 0 more than elementary, 1 elementary
data_raw$education <- if_else(data_raw$education!=0,0,1)

#remove rows in which bmi combined is missing
data_raw <- data_raw[-which(is.na(data_raw$bmi_combined)),]

#rename covariates for the sake of generalization
new_names <- c("bmi_combined" = "cov1", "Smoking" = "cov2", "education"="cov3", "age_wave1"="age")
names(data_raw)[names(data_raw) %in% names(new_names)] <- new_names[names(data_raw)[names(data_raw) %in% names(new_names)]]

data_long <- msprep(data = data_raw, trans = tmat, 
                    time = c(NA, "onset_age", "death_time"), 
                    status = c(NA, "onset", "dead"), 
                    keep = c("age","cov1", "cov2", "cov3"),
                    id="lopnr")

data_long$Tstart[data_long$trans<3] <- data_long$Tstart[data_long$trans<3] +data_long$age[data_long$trans<3]
data_long$time <- data_long$Tstop-data_long$Tstart

# distribution cov 1
plot(data_raw$age, data_raw$cov1)
qqnorm( data_raw$cov1)
qqline(data_raw$cov1)

qqnorm( log(data_raw$cov1))
qqline(log(data_raw$cov1))

meanlog <- mean(log(data_raw$cov1))
sdlog <- sd(log(data_raw$cov1))

#########################
# Fitting a Markov Model
#########################

#fitting a parametric model for each transition adjust for covariates at their baseline level

n_trans <- max(tmat, na.rm = TRUE)
fits_wei <- vector(mode = "list", length = n_trans)
for (i in 1:3){
  fits_wei[[i]] <- flexsurvreg(Surv(Tstart, Tstop, status) ~ cov1+cov2+cov3,  
                               data = subset(data_long, trans == i),
                               dist = "gompertz")
} 

# 1 is dementia-free-> dementia
# 2 is dementia-free-> death
# 3 is dementia-> death

#diagnostic of the cumulative risk for each transition

# plot(fits_wei[[1]],type="cumhaz", xlim=c(60,100))
# plot(fits_wei[[2]],type="cumhaz", xlim=c(60,100))
# plot(fits_wei[[3]],type="cumhaz", xlim=(c(min(data_raw$age_death[data_raw$onset==1], na.rm = T),100)))

#changing coefficients value
# cov 1 continuous such as blood pressure, bmi, high values increases hazards 1->3 and 2->3
# cov 2 binary 15% prevalence such as diabetes or smoking, its presence increases hazards on all transitions
# cov 3 binary 30% prevalence such as low education, its presence increases hazards on 1->2 and 2->3

#we set only factors the increase the risk of transition, no protective factors

fits_wei[[1]]$coefficients[3:5] <- c(0,0.7,0.4) # dementia-free -> dementia cov2 2*risk, cov3 1.5*risk
fits_wei[[2]]$coefficients[3:5] <- c(0.4,0.7,0.2) # dementia-free -> death cov1 1.5*risk, cov2 2*risk, cov3 1.2*risk
fits_wei[[3]]$coefficients[3:5] <- c(0.4,0,0.4) # dementia ->death  cov1 1.5*risk, cov3 1.5*risk

fits_wei[[1]]$res.t[3:5,1] <- c(0,0.7,0.4) 
fits_wei[[2]]$res.t[3:5,1] <- c(0.4,0.7,0.2) 
fits_wei[[3]]$res.t[3:5,1] <- c(0.4,0,0.4) 


fits_wei[[1]]$res[3:5,1] <- c(0,0.7,0.4) 
fits_wei[[2]]$res[3:5,1] <- c(0.4,0.7,0.2) 
fits_wei[[3]]$res[3:5,1] <- c(0.4,0,0.4) 

ground_truth_params <- matrix(0, nrow = 3, ncol = 5)
param_names <- names(fits_wei[[1]]$coefficients)
names(ground_truth_params) <- param_names

for (i in 1:3){
  for (j in 1:5){
    ground_truth_params[i,j] <- round(fits_wei[[i]] $coefficients[j],5)
  }
}

save(ground_truth_params, fits_wei, file = "ground_truthMM.RData")
  

#############################
# Fitting a Semi-Markov Model
#############################

n_trans <- max(tmat, na.rm = TRUE)
fits_wei <- vector(mode = "list", length = n_trans)
for (i in 1:2){
  fits_wei[[i]] <- flexsurvreg(Surv(Tstart, Tstop, status) ~ cov1+cov2+cov3,  
                               data = subset(data_long, trans == i),
                               dist = "gompertz")
}
fits_wei[[3]] <- flexsurvreg(Surv(time, status) ~ cov1+cov2+cov3,  
                               data = subset(data_long, trans == 3),
                               dist = "gompertz")

# for (i in 1:3){
#   print(fits_wei[[i]])
# }

#diagnostic of the cumulative risk for each transition Semi

# plot(fits_wei[[1]],type="cumhaz", xlim=c(60,100))
# plot(fits_wei[[2]],type="cumhaz", xlim=c(60,100))
# plot(fits_wei[[3]],type="cumhaz", xlim=c(0,20))

fits_wei[[1]]$coefficients[3:5] <- c(0,0.7,0.4) # dementia-free -> dementia cov2 2*risk, cov3 1.5*risk
fits_wei[[2]]$coefficients[3:5] <- c(0.4,0.7,0.2) # dementia-free -> death cov1 1.5*risk, cov2 2*risk, cov3 1.2*risk
fits_wei[[3]]$coefficients[3:5] <- c(0.4,0,0.4) # dementia ->death  cov1 1.5*risk, cov3 1.5*risk

fits_wei[[1]]$res.t[3:5,1] <- c(0,0.7,0.4) 
fits_wei[[2]]$res.t[3:5,1] <- c(0.4,0.7,0.2) 
fits_wei[[3]]$res.t[3:5,1] <- c(0.4,0,0.4) 

fits_wei[[1]]$res[3:5,1] <- c(0,0.7,0.4) 
fits_wei[[2]]$res[3:5,1] <- c(0.4,0.7,0.2) 
fits_wei[[3]]$res[3:5,1] <- c(0.4,0,0.4) 


ground_truth_params <- matrix(0, nrow = 3, ncol = 5)
param_names <- names(fits_wei[[1]]$coefficients)
names(ground_truth_params) <- param_names

for (i in 1:3){
  for (j in 1:5){
    ground_truth_params[i,j] <- round(fits_wei[[i]] $coefficients[j],5)
  }
}

save(ground_truth_params, fits_wei, file = "ground_truthSM.RData")

############################
# Data Simulation with hesim
############################

# set model_type "forward" if you want Markov assumption to hold, "mix" otherwise and upload data accordingly

n_pats <- 500
model_type <- "mix"
model <- load("ground_truthSM.RData")

n_sim <- 10
#set.seed(seed)
seeds <- c(555:654)

followup <- 20
dataset <- vector("list", length(n_sim))

# *dataset* this a vector of n_sim elements in which each element is a list of 5 data sets following 
# 5 different observation schemes with fixed n_pats and model_type.
# *probs* this a vector of n_sim elements in which each element represents a dataframe containg the occupancy
# probabilities for that specific simulation


#  paralleling over seeds: 500 patients for 10 seeds we go from 45 to 32 seconds 
#  how to parallelise over patients id to introduce obs scheme?

start_time <- Sys.time()
dataset <- lapply(1:n_sim, function(i) {
  run_simulation(n_pats, fits_wei, model_type, seeds[i], followup, meanlog, sdlog)
})

end_time <- Sys.time()
elapsed_time <- end_time - start_time
elapsed_time

# doesn't really improve the performances
start_time <- Sys.time()
dataset <- mclapply(1:n_sim, function(i){
  run_simulation(n_pats, fits_wei, model_type, seeds[i], followup, meanlog, sdlog)
}, mc.cores= detectCores())
end_time <- Sys.time()
elapsed_time <- end_time - start_time
elapsed_time


setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/Simulation/saved_data")
for (i in 1:n_sim) {
  filename <- paste0("simulation100K_MM_seed_", i, ".RData")
  data <- dataset[[i]]
  save(data, file = filename)
}


########################
# Diagnostic of results
########################

#load("simulation5K_MM_oneseed.RData")
load("simulation5K_MM_10seed.RData")

n_pats <- nrow(dataset[[1]][[1]])

tables <- vector("list", length(dataset))

for (i in 1:length(dataset)){
  list_temp <- dataset[[i]]
  EO <- sum(list_temp[[1]]$onset)
  
  res_1 <- list_temp[[2]] %>%
    group_by(patient_id) %>%              
    summarize(onset_sum = any(onset == 1)) %>%  
    summarize(total_onset = sum(onset_sum, na.rm = TRUE)) %>%
    pull(total_onset)
  res_2 <- list_temp[[3]] %>%
    group_by(patient_id) %>%              
    summarize(onset_sum = any(onset == 1)) %>%  
    summarize(total_onset = sum(onset_sum, na.rm = TRUE)) %>%
    pull(total_onset)
  res_3 <- list_temp[[4]] %>%
    group_by(patient_id) %>%              
    summarize(onset_sum = any(onset == 1)) %>%  
    summarize(total_onset = sum(onset_sum, na.rm = TRUE)) %>%
    pull(total_onset)
  res_4 <- list_temp[[5]] %>%
    group_by(patient_id) %>%              
    summarize(onset_sum = any(onset == 1)) %>%  
    summarize(total_onset = sum(onset_sum, na.rm = TRUE)) %>%
    pull(total_onset)
    
  tables[[i]] <- data.frame(
    dementia_onset = c( EO, res_1, res_2, res_3, res_4),
    row.names = c("Exactly observed", "Scheme 1y", "Scheme 3y", "Scheme Snac-k", "Scheme UkBiobank"))
}

onset_comparison <- do.call(cbind, tables)
onset_comparison_prop <- onset_comparison
onset_comparison_prop[1:5,] <- onset_comparison_prop[1:5,]/n_pats*100
onset_comparison_diff <- apply(onset_comparison_prop, 2, function(x) x[1]-x)
onset_comparison_prop$mean_onset <- rowMeans(onset_comparison_diff)


# onset_comparison[3,]>onset_comparison[2,] # schemes A always spots more onset than scheme B thanks to narrow interval
# onset_comparison[4,]>onset_comparison[3,] # Snac-k doesn't always perform worst than B
# onset_comparison[5,]>onset_comparison[4,] # Ukbiobank doesn't always perform worst than Snack scheme

tables_times <- vector("list", length(dataset))

for (i in 1:length(dataset)){
  list_temp <- dataset[[i]]
  EO <- list_temp[[1]]
  age_onset_1y <- list_temp[[2]] %>%
    group_by(patient_id) %>%
    summarise(onset_age = first(unique(onset_age)))
  age_onset_3y <- list_temp[[3]] %>%
    group_by(patient_id) %>%
    summarise(onset_age = first(unique(onset_age)))
  age_onset_snack <- list_temp[[4]] %>%
    group_by(patient_id) %>%
    summarise(onset_age = first(unique(onset_age)))
  age_onset_ukbiobank <- list_temp[[5]] %>%
    group_by(patient_id) %>%
    summarise(onset_age = first(unique(onset_age)))
  EO$one_y <- age_onset_1y[,2]
  EO$three_y <- age_onset_3y[,2]
  EO$snack <- age_onset_snack[,2]
  EO$ukbiobank <- age_onset_ukbiobank[,2]
  
  diff_age_onset_1y <- EO$one_y-EO$onset_age
  diff_age_onset_3y <- EO$three_y-EO$onset_age
  diff_age_onset_snack <- EO$snack-EO$onset_age
  diff_age_onset_ukbiobank <- EO$ukbiobank-EO$onset_age
  
  tables_times[[i]] <- data.frame(
    mean_delay = c( colMeans(diff_age_onset_1y), colMeans(diff_age_onset_3y), colMeans(diff_age_onset_snack), colMeans(diff_age_onset_ukbiobank)),
    row.names = c("Scheme 1y", "Scheme 3y", "Scheme Snac-k", "Scheme UkBiobank"))
}

onset_delay_comparison <- do.call(cbind, tables_times)
onset_delay_comparison$average_dealy <- rowMeans(onset_delay_comparison)
onset_delay_comparison$average_dealy_days <- ceiling(onset_delay_comparison$average_dealy*12*30)
  
# for (i in 1:10){
#   print(
#     ggplot(probs[[i]], aes(x = t, y = prob)) +
#       geom_line(size = 1) +
#       facet_grid("state_id", scales="free_y")+
#       labs(title = "State occupancy over Time",
#            x = "Age",
#            y = "Probability") +
#       theme_minimal() 
#   )
# }

# salvare data di onset per ciascun paziente in ciascun datset 
# check everything is fine ex patient 24
# data <- dataset [[1]]
# EO <- data[[1]]
# data1 <- data[[2]]
# data3 <- data[[3]]
# data_snack <- data[[4]]
# data_Uk <- data[[5]]

## Compare true parameters to the onse found fitting model over EO dataset
data <- dataset[[1]]
EO_data <- data[[1]] 


check <- EO_data[EO_data$death_time<EO_data$age,] #no weird behaviour
print(max(EO_data$death_time))
table(EO_data$death_time>100)
table(data_raw$death_time>100)

tmat <- mstate::transMat(x = list(c(2, 3),c(3),c()), names = c("Dementia-free","Dementia", "Death")) 


data_long_EO <- msprep(data = EO_data, trans = tmat, 
                       time = c(NA, "onset_age", "death_time"), 
                       status = c(NA, "onset", "dead"), 
                       keep = c("age","cov1", "cov2", "cov3"),
                       id="patient_id")

data_long_EO$Tstart[data_long_EO$trans<3] <- data_long_EO$Tstart[data_long_EO$trans<3] +data_long_EO$age[data_long_EO$trans<3]

data_long_EO$time <- data_long_EO$Tstop-data_long_EO$Tstart

n_trans <- max(tmat, na.rm = TRUE)
fits_wei_EO <- vector(mode = "list", length = n_trans)

for (i in 1:3){
  fits_wei_EO[[i]] <- flexsurvreg(Surv(Tstart, Tstop, status) ~ cov1+cov2+cov3,  
                                  data = subset(data_long_EO, trans == i),
                                  dist = "gompertz")
}
# fits_wei_EO[[3]] <- flexsurvreg(Surv(time, status) ~ cov1+cov2+cov3,  
#                                 data = subset(data_long_EO, trans == 3),
#                                 dist = "gompertz")


params_EO <- matrix(0, nrow = n_trans, ncol = 5)
param_names <- names(fits_wei_EO[[1]]$coefficients)
colnames(params_EO) <- param_names

for (i in 1:3){
  for (j in 1:5){
    params_EO[i,j] <- fits_wei_EO[[i]] $coefficients[j]
  }
}      

params_EO
ground_truth_params

