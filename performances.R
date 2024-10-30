####################################
# Upload library and data
####################################

library(fs)
library(elect)
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

setwd(here())

load("./wrapper_MM/ground_truthMM.RData")
source("./wrapper_MM/functions_performance/compute_bias.R")
source("./wrapper_MM/functions_performance/hazards_mine.R")
source("./wrapper_MM/functions_performance/run_performance_bias.R")
source("./wrapper_MM/functions_performance/run_performance_coverage.R")
source("./wrapper_MM/functions_performance/compute_CI.R")
source("./wrapper_MM/functions_performance/compute_coverage.R")
source("./wrapper_MM/functions_performance/get_params_nhm.R")
source("./wrapper_MM/functions_performance/mean_bias_comparison.R")
source("./wrapper_MM/functions_performance/mean_coverage_comparison.R")
source("./wrapper_MM/functions_performance/check_convergence.R")
source("./wrapper_MM/functions_performance/level_convergence.R")
source("./wrapper_MM/functions_performance/wrapper_convergence.R")
source("./wrapper_MM/functions_performance/plot_convergence.R")
source("./wrapper_MM/functions_performance/plot_bias.R")
source("./wrapper_MM/functions_performance/plot_coverage.R")
source("./wrapper_MM/functions_performance/extract_comp_time.R")
source("./wrapper_MM/functions_performance/plot_boxplot.R")
source("./wrapper_MM/functions_performance/plot_ct.R")
source("./wrapper_MM/functions_performance/prepare_data_boxplot.R")
source("./wrapper_MM/functions_performance/gt_flexsurv.R")

# this code has to be run over each different sample size, is not taken as parameter !
# select number of patients and core to use 

n_pats <- 500
cores <- 4
#scheme <- 4

#######################################################
# fitting parametric model over exactly observed data
#######################################################

# for (seed in 1:100){
#   gt_flexsurv(n_pats, seed)
# }

######################
# check of convergence
######################

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


#######################
# bias comparison
#######################

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
  
  plan(multisession, workers = cores) 
  results_list <- future_lapply(1:100, function(seed) {
    temp_results <- run_performance_bias(n_pats, scheme, seed, combined_cov[[scheme - 1]])
    return(temp_results)
  })
  
  results_bias <- do.call(rbind, results_list)
  
  bias_all_schemes[[scheme-1]] <- results_bias
}

res_bias <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  res_bias[[scheme-1]] <- mean_bias_comparison(bias_all_schemes, scheme)
}



#######################
# coverage comparison
#######################

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


#####################
# life expectancy
#####################


#####################
# computational time
#####################

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



##########
# Plots
##########

titles <- c("Convergence for scheme 1y", "Convergence for scheme 3y", "Convergence for Snac-k", "Convergence for UkBiobank")
plot_convergence(2, titles)
plot_convergence(3, titles)
plot_convergence(4, titles)
plot_convergence(5, titles)

titles <- c("Bias for scheme 1y", "Bias for scheme 3y", "Bias for Snac-k", "Bias for UkBiobank")
plot_bias(2, titles)
plot_bias(3, titles)
plot_bias(4, titles)
plot_bias(5, titles)

# keep in mind that in these estimates of the bias are accounted also those models for which convergence was reached but
# but not to the optimum, so that might increase bias extremely
# yess much better removing them from the computation
# coxph performance for wide intervals performs really bad on some datasets

# choose scheme 
scheme <- 4 #2:5
list_data <- prepare_data_boxplot(scheme)

parameters <- c("cov1","cov2","cov3","rate","shape")
plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[1])
plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[2])
plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[3])
plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[4])
plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[5])


titles <- c("95% Coverage for scheme 1y", "95% Coverage for scheme 3y", "95% Coverage for Snac-k", "95% Coverage for UkBiobank")
plot_coverage(2, titles)
plot_coverage(3, titles)
plot_coverage(4, titles)
plot_coverage(5, titles)

titles <- c("Mean CT for scheme 1y", "Mean CT for scheme 3y", "Mean CT  for Snac-k", "Mean CT  for UkBiobank")
plot_ct(2, titles, combined_cov[[1]])
plot_ct(3, titles, combined_cov[[2]])
plot_ct(4, titles, combined_cov[[3]])
plot_ct(5, titles, combined_cov[[4]])

model_dir <- here()
setwd(model_dir)
if (n_pats==500){
  model_dir <- paste0("wrapper_MM/saved_performance_500")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)
} else if (n_pats==2000){
  model_dir <- paste0("wrapper_MM/saved_performance_2K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)
} else if (n_pats==5000){
  model_dir <- paste0("wrapper_MM/saved_performance_5K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)  
} else if (n_pats==10000){
  model_dir <- paste0("wrapper_MM/saved_performance_10K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)  }

save(combined_cov, file = file.path(model_dir,"convergence.RData"))
save(bias_all_schemes, file = file.path(model_dir,"all_bias.RData"))
save(res_bias, file = file.path(model_dir,"bias.RData"))
save(res_cov, file = file.path(model_dir,"95%coverage.RData"))
save(ct_all_schemes, file = file.path(model_dir,"comp_time.RData"))
 
