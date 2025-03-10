################################
# Upload libraries & functions
################################

library(dplyr)
library(haven)
library(mstate)
library(flexsurv)
library(lubridate)
library(hesim)
library(data.table)
library(ggplot2)
library(future)
library(future.apply)
library(parallel)
library(here())

setwd(here())

source("./simulation_functions/prepare_data.R")
source("./simulation_functions/simulation.R")
source("./simulation_functions/apply_right_censoring.R")
source("./simulation_functions/interval_censoring.R")
source("./simulation_functions/run_simulation.R")
source("./simulation_functions/interval_censoring_id.R")

load("mean_cov1.RData")
load("sd_cov1.RData")


############################
# Data Simulation with hesim
############################

n_pats <- 100000
model_type <- "mix"
model <- load("./wrapper_SM/ground_truthSM.RData")

n_sim <- 100
seeds <- c(555:654)

followup <- 20
dataset_all_SM <- vector("list", length(n_sim))

start_time <- Sys.time()
dataset_all_SM <- mclapply(1:n_sim, function(i){
  run_simulation(n_pats, fits_wei, model_type, seeds[i], followup, meanlog, sdlog)
}, mc.cores= detectCores())
end_time <- Sys.time()
elapsed_time <- end_time - start_time
elapsed_time


dir.create("Simulated_data_SM", showWarnings = FALSE, recursive= T)


#######################################################################
dataset_all_SM_500 <- dataset_all_SM
n_pats = 5

for (i in 1:n_sim){
  for (j in 1:5){
    data_tot <- dataset_all_SM_500[[i]]
    data <- data_tot[[1]]
    tiny_data <- data %>%
      group_by(onset,dead) %>%
      sample_n(size = round(n_pats * n()/nrow(data))) %>%
      ungroup()
    keep_id <- tiny_data$patient_id
    temp <- dataset_all_SM_500[[i]][[j]]
    dataset_all_SM_500[[i]][[j]] <- temp[temp$patient_id %in% keep_id, ]
  }}
save(dataset_all_SM_500, file="Simulated_data_SM_prova/simulation500_SM_all.Rdata")

#######################################################################
dataset_all_SM_5K <- dataset_all_SM
n_pats = 5000

for (i in 1:n_sim){
  for (j in 1:5){
    data_tot <- dataset_all_SM_5K[[i]]
    data <- data_tot[[1]]
    tiny_data <- data %>%
      group_by(onset,dead) %>%
      sample_n(size = round(n_pats * n()/nrow(data))) %>%
      ungroup()
    keep_id <- tiny_data$patient_id
    temp <- dataset_all_SM_5K[[i]][[j]]
    dataset_all_SM_5K[[i]][[j]] <- temp[temp$patient_id %in% keep_id, ]
  }}
save(dataset_all_SM_5K, file="Simulated_data_SM/simulation5K_SM_all.Rdata")

#######################################################################
dataset_all_SM_2K <- dataset_all_SM
n_pats = 2000

for (i in 1:n_sim){
  for (j in 1:5){
    data_tot <- dataset_all_SM_2K[[i]]
    data <- data_tot[[1]]
    tiny_data <- data %>%
      group_by(onset,dead) %>%
      sample_n(size = round(n_pats * n()/nrow(data))) %>%
      ungroup()
    keep_id <- tiny_data$patient_id
    temp <- dataset_all_SM_2K[[i]][[j]]
    dataset_all_SM_2K[[i]][[j]] <- temp[temp$patient_id %in% keep_id, ]
  }}
save(dataset_all_SM_2K, file="Simulated_data_SM/simulation2K_SM_all.Rdata")

#######################################################################

dataset_all_SM_10K <- dataset_all_SM
n_pats = 10000

for (i in 1:n_sim){
  for (j in 1:5){
    data_tot <- dataset_all_SM_10K[[i]]
    data <- data_tot[[1]]
    tiny_data <- data %>%
      group_by(onset,dead) %>%
      sample_n(size = round(n_pats * n()/nrow(data))) %>%
      ungroup()
    keep_id <- tiny_data$patient_id
    temp <- dataset_all_SM_10K[[i]][[j]]
    dataset_all_SM_10K[[i]][[j]] <- temp[temp$patient_id %in% keep_id, ]
  }}
save(dataset_all_SM_10K, file="Simulated_data_SM/simulation10K_SM_all.Rdata")

#######################################################################

