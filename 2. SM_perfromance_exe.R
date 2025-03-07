####################################
# Upload library and data
####################################

library(fs)
library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(here)
library(future)
library(future.apply)
library(mstate)
library(flexsurv)
library(deSolve)
library(hesim)
library(kableExtra)
library(webshot)

setwd(here())

load("./wrapper_SM/ground_truthSM.RData")

source_files <- c(
  "./wrapper_SM/functions_performance/check_convergence.R",
  "./wrapper_SM/functions_performance/compute_bias.R",
  "./wrapper_SM/functions_performance/compute_CI.R",
  "./wrapper_SM/functions_performance/compute_coverage.R",
  "./wrapper_SM/functions_performance/compute_power.R",
  "./wrapper_SM/functions_performance/computing_life_expectancy.R",
  "./wrapper_SM/functions_performance/extract_comp_time.R",
  "./wrapper_SM/functions_performance/get_params_nhm.R",
  "./wrapper_SM/functions_performance/gt_flexsurv.R",
  "./wrapper_SM/functions_performance/ic_comparison.R",
  "./wrapper_SM/functions_performance/level_convergence.R",
  "./wrapper_SM/functions_performance/mean_bias_comparison.R",
  "./wrapper_SM/functions_performance/mean_coverage_comparison.R",
  "./wrapper_SM/functions_performance/mean_lfe_comparison.R",
  "./wrapper_SM/functions_performance/mean_power.R",
  "./wrapper_SM/functions_performance/mean_width_ic.R",
  "./wrapper_SM/functions_performance/power_categorical.R",
  "./wrapper_SM/functions_performance/run_performance_bias.R",
  "./wrapper_SM/functions_performance/run_performance_coverage.R",
  "./wrapper_SM/functions_performance/simulation_probs.R",
  "./wrapper_SM/functions_performance/table_power.R",
  "./wrapper_SM/functions_performance/type_1_error.R",
  "./wrapper_SM/functions_performance/wrapper_convergence.R"
)

lapply(source_files, source)

n_pats <- 10000
cores <- detectCores()


###################
# load quantities #
###################

setwd(here())

if (n_pats==500){
  load("./Simulated_data_SM/simulation500_SM_all.RData")
  data <- dataset_all_SM_500
  path <- "./wrapper_SM/saved_performance_500/"
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  lapply(files, load, .GlobalEnv)
}else if (n_pats==2000){
  load("./Simulated_data_SM/simulation2K_SM_all.RData")
  data <- dataset_all_SM_2K
  path <- "./wrapper_SM/saved_performance_2K/"
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  lapply(files, load, .GlobalEnv)
}else if (n_pats==5000){
  load("./Simulated_data_SM/simulation5K_SM_all.RData")
  data <- dataset_all_SM_5K
  path <- "./wrapper_SM/saved_performance_5K/"
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  lapply(files, load, .GlobalEnv)
}else if (n_pats==10000){
  load("./Simulated_data_SM/simulation10K_SM_all.RData")
  data <- dataset_all_SM_10K
  path <- "./wrapper_SM/saved_performance_10K/"
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  lapply(files, load, .GlobalEnv)
}

# directory to save things
model_dir <- here()
setwd(model_dir)
if (n_pats==500){
  model_dir <- paste0("wrapper_SM/saved_performance_500")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)
} else if (n_pats==2000){
  model_dir <- paste0("wrapper_SM/saved_performance_2K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)
} else if (n_pats==5000){
  model_dir <- paste0("wrapper_SM/saved_performance_5K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)  
} else if (n_pats==10000){
  model_dir <- paste0("wrapper_SM/saved_performance_10K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)  }


#######################################################
# fitting parametric model over exactly observed data #
#######################################################

plan(multisession, workers = cores)
future_lapply(1:100, function(seed) {gt_flexsurv(n_pats, seed)})

###############
# convergence #
###############

temp <- vector(mode = "list", length = 4)
convergence_schemes <- vector(mode = "list", length = 4)
hessian_schemes <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  temp[[scheme-1]] <- wrapper_convergence(n_pats, scheme, seed ) 
  convergence_schemes[[scheme-1]] <- temp[[scheme-1]][[1]]
  hessian_schemes[[scheme-1]] <- temp[[scheme-1]][[2]]
}

# we set combined_conv in the following way, for each method
# 0 if algorithm criteria of convergence were not met 
# 1 if the algorithm converged but not to the optimum, so no hessian exists (a bit different for coxph)
# 2 if the algorithm converged to the optimum

combined_cov <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  combined_cov[[scheme-1]] <- level_convergence(scheme)
}

setwd(here())
save(combined_cov, file = file.path(model_dir,"convergence.RData"))

#######################################
#  absolute bias of covariates effect #
#######################################

bias_all_schemes <- vector(mode = "list", length = 4)
estimates <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  results <- data.frame(
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


  results_list <- future_lapply(1:100, function(seed) {
    temp_results <- run_performance_bias(n_pats, scheme, seed, combined_cov[[scheme - 1]])
    print(seed)
    return(temp_results)
  })
  
  temp_bias <- list()
  temp_est <- list()
  for (seed in 1:100){
    temp_bias[[seed]] <-results_list[[seed]][[1]]
    temp_est[[seed]] <-results_list[[seed]][[2]]
  }
  
  results <- do.call(rbind, temp_bias)
  bias_all_schemes[[scheme-1]] <- results
  
  results<- do.call(rbind, temp_est)
  estimates[[scheme-1]] <- results
}

res_bias <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  res_bias[[scheme-1]] <- mean_bias_comparison(bias_all_schemes, scheme)
}

mean_estimates <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  mean_estimates[[scheme-1]] <- mean_bias_comparison(estimates, scheme)
}

setwd(here())
save(estimates, file = file.path(model_dir,"all_estimates.RData")) 
save(mean_estimates, file = file.path(model_dir,"mean_estimates.RData"))
save(bias_all_schemes, file = file.path(model_dir,"bias_all.RData"))
save(res_bias, file = file.path(model_dir,"res_bias.RData"))


# ######################################
#  95 % coverage of covariates effect  #
# ######################################

coverage_all_schemes <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  results_coverage <- data.frame(
    rate = numeric(0),
    shape = numeric(0),
    cov1 = numeric(0),
    cov2 = numeric(0),
    cov3 = numeric(0),
    model = character(0),
    seed = integer(0)
  )
  
  plan(multisession, workers = cores) 
  results_list <- future_lapply(1:100, function(seed) {
    temp_results <- run_performance_coverage(n_pats, scheme, seed, combined_cov[[scheme-1]])
    return(temp_results)
  })
  
  results_coverage <- do.call(rbind, results_list)
  
  coverage_all_schemes[[scheme-1]] <- results_coverage
  
}

res_cov <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  res_cov[[scheme-1]] <- mean_coverage_comparison(coverage_all_schemes, scheme)
}

res_cov_ic <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  res_cov_ic[[scheme-1]] <- ic_comparison(coverage_all_schemes, scheme)
}

save(coverage_all_schemes, file = file.path(model_dir,"all_coverage.RData"))
save(res_cov_ic, file = file.path(model_dir,"95%coverage.RData"))


###################
# life expectancy #
###################

cores_lfe=4
lfe_bias <- vector(mode = "list", length = 4)
lfe_estimates <- vector(mode = "list", length = 4)
covs <-  data.frame(
  cov1 = 0.35, 
  cov2 = 0.15,
  cov3 = 0.30
)

for (scheme in 2:5){
  results <- data.frame(
    lfe = numeric(0),
    model = character(0),
    seed = integer(0)
  )
  
  plan(multisession, workers = cores) 
  results_list <- future_lapply(1:10, function(seed) {
    t_start <- 60
    temp_results <- computing_life_expectancy(n_pats, scheme, seed, combined_cov[[scheme-1]], t_start, covs)
    return(temp_results)
  })
  
  temp_bias <- list()
  temp_est <- list()
  for (seed in 1:10){
    temp_est[[seed]] <-results_list[[seed]][[1]]
    temp_bias[[seed]] <-results_list[[seed]][[2]]
  }
  
  results <- do.call(rbind, temp_bias)
  lfe_bias[[scheme-1]] <- results
  
  results<- do.call(rbind, temp_est)
  lfe_estimates[[scheme-1]] <- results
  
  print(scheme)
}

res_bias_lfe <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  res_bias_lfe[[scheme-1]] <- mean_lfe_comparison(lfe_bias, scheme)
}

mean_estimates_lfe<- vector(mode = "list", length = 4)
for (scheme in 2:5){
  mean_estimates_lfe[[scheme-1]] <- mean_lfe_comparison(lfe_estimates, scheme)
}

setwd(here())
save(mean_estimates_lfe, file = file.path(model_dir,"mean_estimates_lfe.RData"))
save(res_bias_lfe, file = file.path(model_dir,"bias_lfe.RData"))



##########################
# Power and type 1 error #
##########################

# In the table I will represent P(p_i<=alpha) i.e. percentage of times for which the covs effect is significant
# the yellow values represent when the covariate was significant in the ground truth 
# H0 is variable not significant
# Type 1 error -> refuse H0|H0 true -> significant|no sign
# Type 2 error ->  accept H0| H0 false -> no sign|sign
# Power -> refuse H0|H0 false -> significant|sign

significancy_all <- vector(mode = "list", length = 4)
significancy <- vector(mode = "list", length = 4)
alpha <- 0.05

for (scheme in 2:5){
  plan(multisession, workers = cores) 
  results_list <- future_lapply(1:100, function(seed) {
    temp_results <- compute_power(n_pats, scheme, seed, combined_cov[[scheme - 1]], alpha)
    return(temp_results)
  })
  
  # for (seed in 1:100){
  #   temp_results <- compute_power(n_pats, scheme, seed, combined_cov[[scheme - 1]], alpha)
  #   print(seed)
  # }

  results<- do.call(rbind, results_list)
  significancy_all[[scheme-1]] <- results
  significancy[[scheme-1]] <- mean_power(results, scheme)
  
}

save(significancy, file = file.path(model_dir,"significancy.RData"))
save(significancy_all, file = file.path(model_dir,"significancy_all.RData"))

significant_covs <- data.frame("cov1"= c(0,1,1), "cov2"= c(1,1,0), "cov3"=c(1,1,1), "transition"=c(1,2,3))

# green type one error
# yellow power
# table_power(significant_covs, significancy, scheme=2)
# table_power(significant_covs, significancy, scheme=3)
# t4 <- table_power(significant_covs, significancy, scheme=4)
# table_power(significant_covs, significancy, scheme=5)
# 
# save_kable(t4, file = "powerM.html")


pw2 <- power_categorical(significant_covs, significancy, scheme=2)
pw3 <- power_categorical(significant_covs, significancy, scheme=3)
pw4 <- power_categorical(significant_covs, significancy, scheme=4)
pw5 <- power_categorical(significant_covs, significancy, scheme=5)

setwd(model_dir)
plots_power <- list(pw2,pw3,pw4,pw5)
save(plots_power, file = "power_500.RData" )
setwd(here())

err2 <- type_1_error(significant_covs, significancy, scheme=2)
err3 <- type_1_error(significant_covs, significancy, scheme=3)
err4 <- type_1_error(significant_covs, significancy, scheme=4)
err5 <- type_1_error(significant_covs, significancy, scheme=5)

setwd(model_dir)
plots_errorI <- list(err2,err3,err4,err5)
save(plots_errorI, file = "typeIerr_500.RData" )
setwd(here())


######################
# computational time #
######################

ct_all_schemes <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  plan(multisession, workers = cores) 
  results_list <- future_lapply(1:100, function(seed) {
    temp_results <- extract_comp_time(n_pats, scheme, seed)
    return(temp_results)
  })
  
  results_ct <- do.call(rbind, results_list)
  
  ct_all_schemes[[scheme-1]] <- results_ct
}

setwd(here())
save(ct_all_schemes, file = file.path(model_dir,"comp_time.RData"))


