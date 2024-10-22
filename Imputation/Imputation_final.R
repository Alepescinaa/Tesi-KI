####################################
# Upload library and data
####################################

knitr::opts_chunk$set(echo = TRUE) 
set.seed(2024)

library(mstate)
library(survival)
library(dplyr)
library(survminer)
library(flexsurv)
library(tidyr)
library(ggplot2)

setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI")
load("./Simulated_data_MM/simulation5K_MM_all.RData")

setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI/Imputation")
source("./functions_imp/run_imputation.R")
source("./functions_imp/compute_CI.R")
source("./functions_imp/compare_CI.R")

data <-dataset_all_MM_2K[[1]] 
n_pats <- nrow(data[[1]])


################################################
# Impute dataset of patients observed every year
################################################

scheme_visits <- data[[2]] 
scheme_data <- data[[2]]

original_index <- unique(scheme_data$patient_id)
scheme_data$last_bfo <- numeric()
for (i in original_index){
  index <- which(scheme_data$onset[scheme_data$patient_id==i]==1)[1]
  scheme_data$last_bfo[scheme_data$patient_id==i] <- scheme_data$age[scheme_data$patient_id==i][index-1]
}
for (i in original_index){
  if(any(scheme_data$onset[scheme_data$patient_id==i]==1))
    scheme_data$onset[scheme_data$patient_id==i] <- 1
}

scheme_data <- scheme_data[order(scheme_data$patient_id),]

row_id <- scheme_visits %>%
     group_by(patient_id) %>%  
     summarise(nrows = n())

scheme_visits$patient_id <- rep(1:n_pats, times=as.numeric(row_id$nrows))
scheme_data$patient_id <- rep(1:n_pats, times=as.numeric(row_id$nrows))

scheme_data <- scheme_data %>%
  group_by(patient_id) %>%
  slice(1) %>%
  ungroup()


n<-length(unique(scheme_data$patient_id))  # number of patients
m<-20             # number of imputation
nIter<-200      # number of iterations
eps<-0.005        # criteria of convergence
nx1<-3                # number of variables in onset model
nx2<-4                # number of variables in proportional hazards death model (include DISEASE variable)


# =================================================
# Initial estimate of hazard of developing dementia
# =================================================
# we fit a cox model for the development of dementia setting start as age of the enter in
# the study and onset_age as the middle point between beginning of the study and when onset is observed
# for those developing dementia and censored time if not observed

temp <- scheme_data  
temp$midage <- 0
temp$midage <- (scheme_data$last_bfo + scheme_data$onset_age)/2
temp$onset_age <- ifelse(temp$onset==1,temp$midage,temp$onset_age )
mod_onset <- coxph(Surv(age,onset_age,onset)~cov1+cov2+cov3, data = temp)
onsethaz<-basehaz(mod_onset,center=FALSE) 
onsetpar<-mod_onset$coefficients


# ====================================
# Initial estimate of hazards of dying 
# ====================================
# we fit a cox model for survival to death distinguishing between those having dementia and those dementia-free
# for both we set the right extreme of the interval at censoring time/ death and time start at midage
# for those having dementia and beginning of the study for the latter
# extra covariate: onset, the presence of dementia is likely tho change the outcome

temp$start<-temp$age
temp$start <- ifelse(temp$onset==1,temp$midage, temp$start )
mod_death<- coxph(Surv(start,death_time,dead)~cov1+cov2+cov3+onset, data = temp)
deathhaz<-basehaz(mod_death,center=FALSE)
deathpar<-mod_death$coefficients

# =========================================
# Imputation of onset status and onset age
# =========================================
pb <- txtProgressBar(min = 0, max = nIter, style = 3)
type <- "forward"


results<- run_imputation(scheme_data, scheme_visits, onsethaz, deathhaz, m, nIter, nx1, nx2, eps, type)
avg_parameters <- results[[1]]
all_fits <- results[[2]]
criteria <- results[[3]]

#iniziato alle 15:20
save(results, file="imputated_params5K.RData")

# ======================
# convergence graph
# ======================

criteria_mod1 <- numeric()
criteria_mod2 <- numeric()
criteria_mod3 <- numeric()

for (i in 1:nIter){
  criteria_mod1[i] <- criteria[[i]][1]
  criteria_mod2[i] <- criteria[[i]][2]
  criteria_mod3[i] <- criteria[[i]][3]
}

criteria_mod1 <- as.numeric(criteria_mod1)
criteria_mod2 <- as.numeric(criteria_mod2)
criteria_mod3 <- as.numeric(criteria_mod3)

data_convergence <- data.frame(iteration = 1:nIter,
                   criteria1 = criteria_mod1,
                   criteria2 = criteria_mod2,
                   criteria3 = criteria_mod3)

data_convergence <- data_convergence[-1,] #avoid first iteration

data_convergence$sum <- rowSums(data_convergence[, 2:4])
data_convergence$mean <- rowMeans(data_convergence[, 2:4])
data_convergence$max <- pmax(data_convergence$criteria1, data_convergence$criteria2, data_convergence$criteria3)

data_convergence_long <- pivot_longer(data_convergence, cols = c(sum, mean, max), 
                          names_to = "quantity", 
                          values_to = "value")

ggplot(data_convergence_long, aes(x = iteration, y = value, color = quantity)) +
  geom_line() +
  geom_point() +
  labs(title = "Variation of Metrics over Iterations",
       x = "Iteration Number",
       y = "Value") +
  theme_minimal()

################################################################
# Fitting a parametric model for each transition over EO dataset
################################################################

EO_data <- data[[1]] 

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

if (type == "forward") {
  for (i in 1:3) {
    fits_wei_EO[[i]] <- flexsurvreg(Surv(Tstart, Tstop, status) ~ cov1 + cov2 + cov3, 
                                 data = subset(data_long_EO, trans == i), 
                                 dist = "gompertz")
  }
} else if (type == "mix") { 
  for (i in 1:2) {
    fits_wei_EO[[i]] <- flexsurvreg(Surv(Tstart, Tstop, status) ~ cov1 + cov2 + cov3, 
                                 data = subset(data_long_EO, trans == i), 
                                 dist = "gompertz")
  }
  fits_wei_EO[[3]] <- flexsurvreg(Surv(time, status) ~ cov1 + cov2 + cov3, 
                               data = subset(data_long_EO, trans == 3), 
                               dist = "gompertz")
}

params_EO <- matrix(0, nrow = n_trans, ncol = 5)
param_names <- names(fits_wei_EO[[1]]$coefficients)
colnames(params_EO) <- param_names

for (i in 1:3){
  for (j in 1:5){
    params_EO[i,j] <- fits_wei_EO[[i]] $coefficients[j]
  }
}



###################
# Diagnostic
###################

if (type == "forward"){
  ground_truth <- load("ground_truthMM.RData")
}else if(type == "mix"){
  ground_truth <- load("ground_truthSM.RData")
}


diff_matrix <- avg_parameters-ground_truth_params
compute_norm <- function(v) {
  sqrt(sum(v^2)) 
}

model_norms <- apply(diff_matrix, 1, compute_norm)
tot_norm <- sum(model_norms)



bias_imputed_gt <- abs(avg_parameters-ground_truth_params)
bias_EO_gt <- abs(params_EO-ground_truth_params)
bias_EO_impute <- abs(avg_parameters-params_EO)

bias_imputed_gt
bias_EO_gt
bias_EO_impute

####################
# Rubin's covariance
####################

## Let's compute parameters covariance with Rubin rule
param_matrix <- list()
for (i in 1:m){
  param_matrix[[i]] <- sapply(all_fits[[i]], coefficients) 
}

# mean within variance (mean of covariance matrix of each imputation)
U_list <- list()
for (i in 1:m){
  U_list[[i]] <- lapply(all_fits[[i]], function(fit) vcov(fit))  
}

U_bar <- list()
for (i in 1:3){
  U_bar[[i]] <- matrix(0, nrow = nrow(U_list[[1]][[i]]), ncol = ncol(U_list[[1]][[i]]))
  
  for (j in 1:m){
    U_bar[[i]] <- U_bar[[i]]+U_list[[j]][[i]]
  }
  U_bar[[i]] <- U_bar[[i]]/m
}

# between variance (variance of parameters estimate)
param_matrix_mod1 <- list()
param_matrix_mod2 <- list()
param_matrix_mod3 <- list()
for (i in 1:m){
  param_matrix_mod1[[i]] <- param_matrix[[i]][,1]
  param_matrix_mod2[[i]] <- param_matrix[[i]][,2]
  param_matrix_mod3[[i]] <- param_matrix[[i]][,3]
}

B1 <- rep(0, 5, 5)
for (j in 1:m) {
  B1 <- B1 + (param_matrix_mod1[[j]] - avg_parameters[1,])%*% t(param_matrix_mod1[[j]] - avg_parameters[1,])
}
B1 <- B1 / (m - 1)

B2 <- rep(0, 5, 5)
for (j in 1:m) {
  B2 <- B2 + (param_matrix_mod1[[j]] - avg_parameters[2,])%*% t(param_matrix_mod1[[j]] - avg_parameters[2,])
}
B2 <- B2 / (m - 1)

B3 <- rep(0, 5 ,5)
for (j in 1:m) {
  B3 <- B3 + (param_matrix_mod1[[j]] - avg_parameters[3,])%*% t(param_matrix_mod1[[j]] - avg_parameters[3,])
}
B3 <- B3 / (m - 1)


#Rubin's rule
T1 <- U_bar[[1]] + (1 + 1/m) * B1
T2 <- U_bar[[2]] + (1 + 1/m) * B2
T3 <- U_bar[[3]] + (1 + 1/m) * B3

# confidence interval computation

CI_mod1<- compute_CI(T1,U_bar[[1]],B1,m,avg_parameters[1,])
CI_mod2<- compute_CI(T2,U_bar[[2]],B2,m,avg_parameters[2,])
CI_mod3<- compute_CI(T3,U_bar[[3]],B3,m,avg_parameters[3,])


########################################
# Intersection CI from EO and imputation
########################################

CI_EO1 <- fits_wei_EO[[1]]$res.t[,2:3]
CI_EO2 <- fits_wei_EO[[2]]$res.t[,2:3]
CI_EO3 <- fits_wei_EO[[3]]$res.t[,2:3]


compare_CI(CI_mod1,CI_EO1)
compare_CI(CI_mod2,CI_EO2)
compare_CI(CI_mod3,CI_EO3)
