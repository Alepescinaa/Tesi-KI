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
library(smms)
library(igraph)

# choose the sample size and upload accordingly the datset, either 500, 2K, 5K
setwd(here())
n_pats <- 500 # CHANGE HERE

source("./wrapper_SM/functions_wrapper/prepare_flex.R")
source("./wrapper_SM/functions_wrapper/prepare_coxph.R")
source("./wrapper_SM/functions_wrapper/prepare_msm.R")
source("./wrapper_SM/functions_wrapper/prepare_imputation.R")
source("./wrapper_SM/functions_wrapper/fit_model.R")
source("./wrapper_SM/functions_wrapper/averaging_params.R")
source("./wrapper_SM/functions_wrapper/run_imputation.R")
source("./wrapper_SM/functions_wrapper/wrapper_functions_SM.R")
source("./wrapper_SM/functions_wrapper/smms_mine.R")


cores <- 4
cores_nhm <- 4

# Windows
if (n_pats == 500){
  load("./Simulated_data_SM/simulation500_SM_all.RData")
  plan(multisession, workers = cores)  
  for (scheme in 2:5) {
    future_lapply(1:100, function(seed) {
      data <- dataset_all_SM_500[[seed]][[scheme]] 
      wrapper_functions_SM(data, n_pats, seed, cores_nhm)
      return(paste("Completed seed:", seed))
    })
  }
} else if (n_pats == 2000){
  load("./Simulated_data_SM/simulation2K_SM_all.RData")
  plan(multisession, workers = cores)  
  for (scheme in 2:5) {
    future_lapply(1:100, function(seed) {
      data <- dataset_all_SM_2K[[seed]][[scheme]] 
      wrapper_functions_SM(data, n_pats, seed, cores_nhm)
      return(paste("Completed seed:", seed))
    })}
}else if (n_pats == 5000){
  load("./Simulated_data_SM/simulation5K_SM_all.RData")
  plan(multisession, workers = cores)  
  for (scheme in 2:5) {
    future_lapply(1:100, function(seed) {
      data <- dataset_all_SM_5K[[seed]][[scheme]] 
      wrapper_functions_SM(data, n_pats, seed, cores_nhm)
      return(paste("Completed seed:", seed))
    })
  }}else if (n_pats == 10000){
    load("./Simulated_data_SM/simulation10K_SM_all.RData")
    plan(multisession, workers = cores)  
    for (scheme in 2:5) {
      future_lapply(1:100, function(seed) {
        data <- dataset_all_SM_10K[[seed]][[scheme]] 
        wrapper_functions_SM(data, n_pats, seed, cores_nhm)
        return(paste("Completed seed:", seed))
      })
    }}


