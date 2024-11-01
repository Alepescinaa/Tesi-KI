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
library(parallel)
library(future)
library(future.apply)
library(here)

# choose the sample size and upload accordingly the datset, either 500, 2K, 5K
setwd(here())
n_pats <- 500 # CHANGE HERE

source("./wrapper_MM/functions_wrapper/prepare_coxph_flex.R")
source("./wrapper_MM/functions_wrapper/prepare_msm.R")
source("./wrapper_MM/functions_wrapper/prepare_imputation.R")
source("./wrapper_MM/functions_wrapper/fit_model.R")
source("./wrapper_MM/functions_wrapper/averaging_params.R")
source("./wrapper_MM/functions_wrapper/run_imputation.R")
source("./wrapper_MM/functions_wrapper/wrapper_functions_MM.R")


cores <- 4
cores_nhm <- 1

# Windows
if (n_pats == 500){
  load("./Simulated_data_MM/simulation500_MM_all.RData")
  plan(multisession, workers = cores)  
  for (scheme in 2:5) {
    future_lapply(1:100, function(seed) {
      data <- dataset_all_MM_500[[seed]][[scheme]] 
      n_pats <- length(unique(data$patient_id))
      wrapper_functions_MM(data, n_pats, seed, cores_nhm)
      return(paste("Completed seed:", seed))
    })
  }
} else if (n_pats == 2000){
  load("./Simulated_data_MM/simulation2K_MM_all.RData")
  plan(multisession, workers = cores)  
  for (scheme in 2:5) {
    future_lapply(1:100, function(seed) {
      data <- dataset_all_MM_2K[[seed]][[scheme]] 
      n_pats <- length(unique(data$patient_id))
      wrapper_functions_MM(data, n_pats, seed, cores_nhm)
      return(paste("Completed seed:", seed))
    })}
  }else if (n_pats == 5000){
    load("./Simulated_data_MM/simulation5K_MM_all.RData")
    plan(multisession, workers = cores)  
    for (scheme in 2:5) {
      future_lapply(1:100, function(seed) {
        data <- dataset_all_MM_5K[[seed]][[scheme]] 
        n_pats <- length(unique(data$patient_id))
        wrapper_functions_MM(data, n_pats, seed, cores_nhm)
        return(paste("Completed seed:", seed))
      })
    }}else if (n_pats == 10000){
      load("./Simulated_data_MM/simulation10K_MM_all.RData")
  plan(multisession, workers = cores)  
  for (scheme in 2:5) {
    future_lapply(1:100, function(seed) {
      data <- dataset_all_MM_10K[[seed]][[scheme]] 
      n_pats <- length(unique(data$patient_id))
      wrapper_functions_MM(data, n_pats, seed, cores_nhm)
      return(paste("Completed seed:", seed))
    })
  }}
  
  
  
