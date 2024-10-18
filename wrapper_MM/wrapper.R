####################################
# Upload library and data
####################################

library(mstate)
library(survival)
library(dplyr)
library(msm)
library(flexsurv)
library(nhm)
library(ggplot2)

setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI")
load("./Simulated_data_MM/simulation500_MM_all.RData")
load("./wrapper_MM/ground_truthMM.RData")

setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI/wrapper_MM")
source("./functions_wrapper/prepare_coxph_flex.R")
source("./functions_wrapper/prepare_msm.R")
source("./functions_wrapper/prepare_imputation.R")
source("./functions_wrapper/run_imputation.R")
source("./functions_performance/compute_bias.R")


seed <- 5 # 1:100
scheme <- 2 # 2:5 (1y,3y,snac-k,ukbiobank)

data <-dataset_all_MM_500[[seed]][[scheme]]
n_pats <- length(unique(data$patient_id))
comp_time <- numeric()

if (scheme==2){
  model_dir <- paste0("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI/wrapper_MM/saved_models_scheme2/seed_", seed)
  dir.create(model_dir, showWarnings = FALSE)
} else if (scheme==3){
  model_dir <- paste0("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI/wrapper_MM/saved_models_scheme3/seed_", seed)
  dir.create(model_dir, showWarnings = FALSE)
} else if (scheme==4){
  model_dir <- paste0("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI/wrapper_MM/saved_models_scheme4/seed_", seed)
  dir.create(model_dir, showWarnings = FALSE)
} else if (scheme==5){
  model_dir <- paste0("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI/wrapper_MM/saved_models_scheme5/seed_", seed)
  dir.create(model_dir, showWarnings = FALSE)
}

#####################
# coxph model
#####################

temp <-  prepare_coxph_flex(data, n_pats)
model_cox <- vector(mode = "list", length = 3)

time_cox<- system.time({
for (i in 1:3) {
  model_cox[[i]] <- coxph(Surv(Tstart,Tstop,status) ~ cov1 + cov2 + cov3, data = subset(temp, trans == i))}
})[3]

comp_time[1] <- as.numeric(round(time_cox,3))

save(model_cox, file = file.path(model_dir,"cox_model.RData"))

######################
# flexsurv model
######################

temp <- prepare_coxph_flex(data, n_pats)
fits_gompertz <- vector(mode = "list", length = 3)

time_gomp<- system.time({
for (i in 1:3) {
  fits_gompertz[[i]] <- flexsurvreg(Surv(Tstart, Tstop, status) ~ cov1 + cov2 + cov3, 
                                  data = subset(temp, trans == i), 
                                  dist = "gompertz")}
  })[3]

comp_time[2] <- as.numeric(round(time_gomp,3))

save(fits_gompertz, file = file.path(model_dir, "flexsurv_model.RData"))

######################
# msm model
######################

temp <- prepare_msm(data) 

Q <- rbind(c(0, 1, 1), 
           c(0, 0, 1),
           c(0, 0, 0)) 

time_msm<- system.time({
model.msm <- msm(state ~ age,
                 subject = patient_id,
                 data = temp, 
                 qmatrix = Q,
                 covariates = ~ cov1 + cov2 + cov3 ,
                 gen.inits = TRUE,
                 control = list(fnscale = 1000, maxit = 500),
                 deathexact = TRUE)
})[3]

comp_time[3] <- as.numeric(round(time_msm,3))

#initial_guess <- qmatrix.msm(model.msm)$estimates
msm_estimates <- model.msm$estimates.t

save(model.msm, file = file.path(model_dir, "msm_model.RData"))

#####################
# msm + age model
#####################

temp <- prepare_msm(data)

Q <- rbind(c(0, 1, 1), 
           c(0, 0, 1),
           c(0, 0, 0)) 

intial_guess <- qmatrix.msm(model.msm)$estimates

time_msm_age<- system.time({
  model.msm_age <- msm(state ~ age, 
                       subject = patient_id,
                       data = temp, 
                       qmatrix = intial_guess,
                       covariates = ~ cov1 + cov2 + cov3 + age ,
                       gen.inits= F,
                       control = list(fnscale = 1000, maxit = 500),
                       center= TRUE,
                       deathexact = TRUE) 
})[3]

# + veloce scalato e con bias leggermente minore scalando dati, inizalizzazione Q velocizza un po'
# system.time({
# min_age <- min(temp$age)
# temp$age <- temp$age - min_age
# 
# model_age.msm_scaled <- msm(state ~ age,
#                      subject = patient_id,
#                      data = temp, 
#                      qmatrix = Q,
#                      covariates = ~ cov1 + cov2 + cov3 + age , #elect wants factors as dummy
#                      gen.inits=TRUE,
#                      control = list(fnscale = 1000),
#                      #censor = 99,
#                      center= FALSE,# requested by elect
#                      deathexact = TRUE) 
# })[3]


# params_msm_age <- matrix(model.msm_age$estimates[4:12], nrow = 3, ncol = 3) #1:3 rate 12:15 age
# params_msm_age <- cbind(params_msm_age, exp(params_msm_age[,1]), exp(params_msm_age[,2]), exp(params_msm_age[,3]))
# colnames(params_msm_age) <- colnames(ground_truth_params)[3:8]
# 
# bias_msm_age <- compute_bias(params_msm_age, ground_truth_params)
# 
# params_msm_age <- matrix(model_age.msm_scaled$estimates[4:12], nrow = 3, ncol = 3) #1:3 rate 12:15 age
# params_msm_age <- cbind(params_msm_age, exp(params_msm_age[,1]), exp(params_msm_age[,2]), exp(params_msm_age[,3]))
# colnames(params_msm_age) <- colnames(ground_truth_params)[3:8]
# 
# bias_msm_age_scaled <- compute_bias(params_msm_age, ground_truth_params)
# 
# bias_msm_age<bias_msm_age_scaled

comp_time[4] <- as.numeric(round(time_msm_age,3))

save(model.msm, file = file.path(model_dir, "model_msm_age.RData"))

######################
# nhm model
######################

temp <- prepare_msm(data)

tmat_1 <- rbind(c(0,1,2),c(0,0,3),rep(0,3))
tmat_2 <- rbind(c(0,4,5),c(0,0,6),rep(0,3))
tmat_3 <- rbind(c(0,7,8),c(0,0,9),rep(0,3))

temp$patient_id <- as.factor(temp$patient_id) 
temp=as.data.frame(temp)

#initial_guess <-  append(msm_estimates, rep(0.5, 3), after = 3)
# we have estimates for rate and covs, so we add initial estimate to 0.5 of shape -> not helping
 
# we set split t 50% percentile
find_splits <- function(age) {
  quantiles <- quantile(age, probs = seq(0, 1, 0.1))
  return(quantiles[-c(1, length(quantiles))])  
}

split_points <- find_splits(temp$age)[6:9]

time_nhm <- system.time({
  object_nhm <- model.nhm(state ~ age,
                          subject = patient_id,
                          data = temp, 
                          trans = tmat_1,
                          nonh = tmat_1,
                          type = "gompertz",
                          covariates = c("cov1", "cov2", "cov3"),
                          covm = list(cov1= tmat_1, cov2=tmat_2, cov3=tmat_3),
                          death = T, 
                          death.states = c(3))
  
  
  model_nhm <- nhm(object_nhm, 
                   gen_inits = TRUE,
                   #initial = initial_guess,
                   score_test = FALSE, 
                   control = nhm.control(ncores = 4, obsinfo = FALSE, coarsen = T, coarsen.vars = c(1), coarsen.lv = 5, splits = split_points))
})[3]

params_nhm <- matrix(model_nhm$par, nrow = 3, ncol = 5) 
params_nhm <- cbind(params_nhm, exp(params_nhm[,3]), exp(params_nhm[,4]), exp(params_nhm[,5]))
colnames(params_nhm) <- colnames(ground_truth_params)

bias_nhm_c5 <- compute_bias(params_nhm, ground_truth_params)

time_nhm_cont <- system.time({
  object_nhm <- model.nhm(state ~ age,
                          subject = patient_id,
                          data = temp, 
                          trans = tmat_1,
                          nonh = tmat_1,
                          type = "gompertz",
                          covariates = c("cov1", "cov2", "cov3"),
                          covm = list(cov1= tmat_1, cov2=tmat_2, cov3=tmat_3),
                          death = T, 
                          death.states = c(3))
  
  
  model_nhm <- nhm(object_nhm, 
                   gen_inits = TRUE,
                   #initial = initial_guess,
                   score_test = FALSE, 
                   control = nhm.control(ncores = 4, obsinfo = FALSE, splits = split_points))
})[3]

params_nhm <- matrix(model_nhm$par, nrow = 3, ncol = 5) 
params_nhm <- cbind(params_nhm, exp(params_nhm[,3]), exp(params_nhm[,4]), exp(params_nhm[,5]))
colnames(params_nhm) <- colnames(ground_truth_params)

bias_nhm_cont <- compute_bias(params_nhm, ground_truth_params)

bias_nhm_cont<bias_nhm_c5
bias_increase_cors <- (bias_nhm_c5-bias_nhm_cont)/bias_nhm_cont*100

save(bias_increase_cors, file="bias_increase_cors.RData")

#both model have automatic initialization and splits
#comp.time 500 pats, no coarsening 2052.836 more than 4 times slower than the version with coarsening
#comp.time 500 pats, coarsening to 5 values 431 errore percentual salvato in bias_increase_cors

comp_time[5] <- as.numeric(round(time_nhm,3))

save(model.msm, file = file.path(model_dir, "model_nhm.RData"))

####################
# imputation model
####################

temp <- prepare_imputation(data)
scheme_data <- temp[[1]]
scheme_visits <- temp[[2]]
m <- 30
type <- "forward"

time_imp <- system.time({
  results_imp<- run_imputation(scheme_data, scheme_visits, m, type)
  avg_parameters <- results_imp[[1]]
  all_fits <- results_imp[[2]]
})[3]

comp_time[6] <- as.numeric(round(time_imp,3))

save(model.msm, file = file.path(model_dir, "results_imp.RData"))

save(comp_time, file = file.path(model_dir, "computational_time.RData"))

