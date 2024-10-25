####################################
# Upload library and data
####################################

library(fs)
library(elect)
library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)

load("ground_truthMM.RData")
source("./functions_performance/compute_bias.R")
source("./functions_performance/hazards_mine.R")
source("./functions_performance/run_performance_bias.R")
source("./functions_performance/run_performance_coverage.R")
source("./functions_performance/compute_CI.R")
source("./functions_performance/compute_coverage.R")
source("./functions_performance/get_params_nhm.R")
source("./functions_performance/mean_bias_comparison.R")
source("./functions_performance/mean_coverage_comparison.R")
source("./functions_performance/check_convergence.R")
source("./functions_performance/wrapper_convergence.R")
source("./functions_performance/plot_convergence.R")

setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI/wrapper_MM")

# this code has to be run over each different sample size, is not taken as parameter !
# select number of patients and core to use 

n_pats <- 500
scheme <-  2
cores <- 4

######################
# check of convergence
######################

convergence_schemes <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  convergence_schemes[[scheme-1]] <- wrapper_convergence(n_pats, scheme, seed )
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
  
  results_list <- mclapply(1:100, function(seed) {
    setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI/wrapper_MM")
    temp_results <- run_performance_bias(n_pats, scheme, seed)
    return(temp_results)  
  }, mc.cores = cores)
  
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
  
  results_list <- mclapply(1:100, function(seed) {
    setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI/wrapper_MM")
    temp_results <- run_performance_coverage(n_pats, scheme, seed)
    return(temp_results)  
  }, mc.cores = cores)
  
  results_coverage <- do.call(rbind, results_list)
  
  coverage_all_schemes[[scheme-1]] <- results_coverage
  
}

res_cov <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  res_cov[[scheme-1]] <- mean_coverage_comparison(coverage_all_schemes, scheme)
}


##########
# Plots
##########

titles <- c("Convergence for scheme 1y", "Convergence for scheme 3y", "Convergence for Snac-k", "Convergence for UkBiobank")

for (scheme in 2:5){
  plot_convergence(scheme, titles)
}
