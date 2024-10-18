####################################
# Upload library and data
####################################

library(fs)

main_dir <- "/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-Ki/wrapper_MM/saved_models_scheme"

scheme <-  2
seed <- 1

if (scheme == 2){
  scheme_dir <- "/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-Ki/wrapper_MM/saved_models_scheme2"
} else if (scheme == 3){
  scheme_dir <- "/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-Ki/wrapper_MM/saved_models_scheme3"
} else if (scheme == 4){
  scheme_dir <- "/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-Ki/wrapper_MM/saved_models_scheme4"
} else if (scheme == 5){
  scheme_dir <- "/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-Ki/wrapper_MM/saved_models_scheme5"
}
  

seed_dir <- file.path(scheme_dir, paste0("seed_", seed))

if (dir.exists(seed_dir)) {
  setwd(seed_dir)
  
  files_to_load <- c("cox_model.RData", 
                     "flexsurv_model.RData", 
                     "msm_model.RData", 
                     "model_msm_age.RData", 
                     "model_nhm.RData", 
                     "results_imp.RData", 
                     "computational_time.RData")
  for (file in files_to_load) {
    if (file.exists(file)) {
      load(file)
    } else {
      warning(paste("File does not exist:", file))
    }
  }
} else {
  warning(paste("Seed directory does not exist:", seed_dir))
}

setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-Ki/wrapper_MM")
load("ground_truthMM.RData")
source("./functions_performance/compute_bias.R")

# =========
# coxph
# =========

params_coxph <- matrix(0, nrow = 3, ncol = 3)
param_names <- names(model_cox[[1]]$coefficients)
colnames(params_coxph) <- param_names

for (i in 1:3){
  for (j in 1:3){
    params_coxph[i,j] <- model_cox[[i]] $coefficients[j]
  }
}
params_coxph <- cbind(params_coxph, exp(params_coxph[,1]), exp(params_coxph[,2]), exp(params_coxph[,3]))
colnames(params_coxph)[4:6] <- c("exp(cov1)", "exp(cov2)", "exp(cov3)")

bias_coxph <- compute_bias(params_coxph, ground_truth_params)

# ============
# flexsurv
# ============

params_flexsurv <- matrix(0, nrow = 3, ncol = 5)
param_names <- names(fits_gompertz[[1]]$coefficients)
colnames(params_flexsurv) <- param_names

for (i in 1:3){
  for (j in 1:5){
    params_flexsurv[i,j] <- fits_gompertz[[i]] $coefficients[j]
  }
}

params_flexsurv <- params_flexsurv[, c(2, 1, 3, 4, 5)]
params_flexsurv <- cbind(params_flexsurv, exp(params_flexsurv[,3]), exp(params_flexsurv[,4]), exp(params_flexsurv[,4]))
colnames(params_flexsurv)[6:8] <- c("exp(cov1)", "exp(cov2)", "exp(cov3)")

bias_flexsurv <- compute_bias(params_flexsurv, ground_truth_params)

# ============
# msm
# ============

params_msm <- matrix(model.msm$estimates[4:12], nrow = 3, ncol = 3) #1:3 rate
params_msm <- cbind(params_msm, exp(params_msm[,1]), exp(params_msm[,2]), exp(params_msm[,3]))
colnames(params_msm) <- colnames(ground_truth_params)[3:8]

bias_msm <- compute_bias(params_msm, ground_truth_params)

# ============
# msm + age
# ============

params_msm_age <- matrix(model.msm_age$estimates[4:12], nrow = 3, ncol = 3) #1:3 rate 12:15 age
params_msm_age <- cbind(params_msm_age, exp(params_msm_age[,1]), exp(params_msm_age[,2]), exp(params_msm_age[,3]))
colnames(params_msm_age) <- colnames(ground_truth_params)[3:8]

bias_msm_age <- compute_bias(params_msm_age, ground_truth_params)

# ============
# nhm
# ============

# ============
# imputation
# ============

params_imp <- results_imp[[1]]
params_imp <- params_imp[, c(2, 1, 3, 4, 5)]
params_imp <- cbind(params_imp, exp(params_imp[,3]), exp(params_imp[,4]), exp(params_imp[,4]))
colnames(params_imp)[6:8] <- c("exp(cov1)", "exp(cov2)", "exp(cov3)")

bias_imputation <- compute_bias(params_imp, ground_truth_params)


