####################################
# Upload library and data
####################################

library(fs)
library(elect)
library(dplyr)

load("ground_truthMM.RData")
source("./functions_performance/compute_bias.R")
source("./functions_performance/hazards_mine.R")
source("./functions_performance/run_performance_bias.R")

n_pats <- 500
scheme <-  2
seed <- 3
cores <- 4

bias_all_schemes <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  results_bias <- data.frame(
    rate = numeric(0),
    shape = numeric(0),
    cov1 = numeric(0),
    cov2 = numeric(0),
    cov3 = numeric(0),
    `exp(cov1)` = numeric(0),
    `exp(cov2)` = numeric(0),
    `exp(cov3)` = numeric(0),
    model = character(0),
    seed = integer(0)
  )
  
  results_list <- mclapply(1:3, function(seed) {
    temp_results <- run_performance(n_pats, scheme, seed)
    return(temp_results)  
  }, mc.cores = cores)
  
  results_bias <- do.call(rbind, results_list)
  
  bias_all_schemes[[scheme-1]] <- results_bias
}

temp <- bias_all_schemes[[1]]
temp <- as.data.frame(temp)
temp <- temp %>%
  mutate(across(1:8, as.numeric))

mean_results <- temp %>%
  group_by(model, transition) %>%
  summarise(
    across(c(rate, shape, cov1, cov2, cov3, `exp(cov1)`, `exp(cov2)`, `exp(cov3)`), 
           ~ mean(.x, na.rm = TRUE)), 
    .groups = 'drop'
  )

# =============
# compare bias
# =============

# colMeans(bias_coxph)
# colMeans(bias_flexsurv)
# colMeans(bias_msm)
# colMeans(bias_msm_age)
# colMeans(bias_nhm)
# colMeans(bias_imputation)
# 
# #imputation performs better then other parametric methods using gompertz 
# colMeans(bias_imputation)<colMeans(bias_flexsurv)
# colMeans(bias_imputation)<colMeans(bias_nhm)
# 
# #introducing age as covariate improves estimates for covariate effect
# colMeans(bias_msm_age)<colMeans(bias_msm)
# colMeans(bias_msm)<colMeans(bias_imputation)
# #nhm performs really poorly

# if (scheme==2){
#   model_dir <- paste0("bias_500/scheme2")
#   dir.create(model_dir, showWarnings = FALSE, recursive= T)
# } else if (scheme==3){
#   model_dir <- paste0("bias_500/scheme3")
#   dir.create(model_dir, showWarnings = FALSE, recursive= T)
# } else if (scheme==4){
#   model_dir <- paste0("bias_500/scheme4")
#   dir.create(model_dir, showWarnings = FALSE, recursive= T)  
# } else if (scheme==5){
#   model_dir <- paste0("bias_500/saved_models_scheme5")
#   dir.create(model_dir, showWarnings = FALSE, recursive= T)  
# }
# 
# setwd(model_dir)
# save(results_bias, file = file.path(model_dir,"results_bias.RData"))
