####################################
# Upload library and data
####################################

library(fs)
library(elect)
library(dplyr)
library(parallel)

load("ground_truthMM.RData")
source("./functions_performance/compute_bias.R")
source("./functions_performance/hazards_mine.R")
source("./functions_performance/run_performance_bias.R")
source("./functions_performance/run_performance_coverage_copia.R")
source("./functions_performance/compute_CI.R")
source("./functions_performance/compute_coverage.R")
source("./functions_performance/get_params_nhm.R")
source("./functions_performance/mean_bias_comparison.R")

setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI/wrapper_MM")

# this code has to be run over each different sample size, is not taken as parameter !
# select number of patients and core to use 

n_pats <- 500
scheme <-  2
cores <- 4


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

res <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  res[[scheme-1]] <- mean_bias_comparison(bias_all_schemes, scheme)
}

#######################
# coverage comparison
#######################

coverage_all_schemes <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  results_bias <- data.frame(
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
# this would be a function comparing mean coverage berween models for each scheme 
temp <- coverage_all_schemes[[1]]
temp <- as.data.frame(temp)
temp <- temp %>%
  mutate(across(1:5, as.numeric))

mean_coverage <- temp %>%
  group_by(model, transition) %>%
  summarise(
    across(c(rate, shape, cov1, cov2, cov3), 
           ~ round(mean(.x, na.rm = TRUE), 3)), 
    .groups = 'drop'
  )



