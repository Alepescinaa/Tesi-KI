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

# choose the sample size and upload accordingly the datset, either 500, 2K, 5K
#setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI")
load("./Simulated_data_MM/simulation500_MM_all.RData")

setwd("./wrapper_MM")
source("./functions_wrapper/prepare_coxph_flex.R")
source("./functions_wrapper/prepare_msm.R")
source("./functions_wrapper/prepare_imputation.R")
source("./functions_wrapper/run_imputation.R")
source("./functions_wrapper/wrapper_functions_MM.R")


# loop
for (scheme in 2:5){
  for (seed in 1:2){
    data <-dataset_all_MM_500[[seed]][[scheme]]
    n_pats <- length(unique(data$patient_id))
    wrapper_functions_MM(data, n_pats, seed, cores_nhm)
    print("models completed for seed_", seed)
  }
}

cores <- 4
cores_nhm <- 1
scheme <- 2

# Mac
for (scheme in 2:5) {
  mclapply(1:2, function(seed) {
    data <- dataset_all_MM_500[[seed]][[scheme]]
    n_pats <- length(unique(data$patient_id))
    wrapper_functions_MM(data, n_pats, seed, cores_nhm)
    }, mc.cores = cores)
}


# Windows
for (scheme in 2:5) {
  cl <- makeCluster(cores)
  clusterExport(cl, varlist = c("dataset_all_MM_500", "wrapper_functions_MM", "scheme", "prepare_coxph_flex", "prepare_msm", "prepare_imputation", "run_imputation" ))
  
  parLapply(cl, 1:2, function(seed) {
    data <- dataset_all_MM_500[[seed]][[scheme]]
    n_pats <- length(unique(data$patient_id))
    wrapper_functions_MM(data, n_pats, seed, cores_nhm)
  })
  
  stopCluster(cl)

}


plan(multisession)  
for (scheme in 2:5) {
  future_lapply(1:2, function(seed) {
    data <- dataset_all_MM_500[[seed]][[scheme]]
    n_pats <- length(unique(data$patient_id))
    wrapper_functions_MM(data, n_pats, seed, cores_nhm)
    #return(paste("Completed seed:", seed))
  })
}
 
