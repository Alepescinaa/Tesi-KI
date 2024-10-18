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

# choose the sample size and upload accordingly the datset, either 500, 2K, 5K
setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI")
load("./Simulated_data_MM/simulation500_MM_all.RData")

setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI/wrapper_MM")
source("./functions_wrapper/prepare_coxph_flex.R")
source("./functions_wrapper/prepare_msm.R")
source("./functions_wrapper/prepare_imputation.R")
source("./functions_wrapper/run_imputation.R")
source("./functions_wrapper/wrapper_functions_MM.R")



for (scheme in 2:5){
  for (seed in 1:3){
    data <-dataset_all_MM_500[[seed]][[scheme]]
    n_pats <- length(unique(data$patient_id))
    wrapper_functions_MM(data,n_pats)
    print("models completed for seed_", seed)
  }
}


# for (scheme in 2:5) {
#   mclapply(1:3, function(seed) {
#     data <- dataset_all_MM_500[[seed]][[scheme]]
#     n_pats <- length(unique(data$patient_id))
#     wrapper_functions_MM(data, n_pats)
#   }, mc.cores = 4) 
# }
#   

