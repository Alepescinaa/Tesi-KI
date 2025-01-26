################################
# Upload libraries & functions
################################

knitr::opts_chunk$set(echo = TRUE)

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

source("./functions/prepare_data.R")
source("./functions/simulation.R")
source("./functions/apply_right_censoring.R")
source("./functions/interval_censoring.R")
source("./functions/run_simulation.R")
source("./functions/interval_censoring_id.R")

################################################
# Upload data and preprocess them in long format
################################################

load("SNAC_data_Alessandra.RData")
data_raw <- SNAC_K_dataset

data_raw <- prepare_data(data_raw) 

tmat <- mstate::transMat(x = list(c(2, 3),c(3),c()), names = c("Dementia-free","Dementia", "Death")) 
data_raw$Smoking <- if_else(data_raw$Smoking==0|data_raw$Smoking==1,0,1) 
data_raw$education <- if_else(data_raw$education!=0,0,1)
data_raw <- data_raw[-which(is.na(data_raw$bmi_combined)),]
new_names <- c("bmi_combined" = "cov1", "Smoking" = "cov2", "education"="cov3", "age_wave1"="age")
names(data_raw)[names(data_raw) %in% names(new_names)] <- new_names[names(data_raw)[names(data_raw) %in% names(new_names)]]

data_long <- msprep(data = data_raw, trans = tmat, 
                    time = c(NA, "onset_age", "death_time"), 
                    status = c(NA, "onset", "dead"), 
                    keep = c("age","cov1", "cov2", "cov3"),
                    id="lopnr")

data_long$Tstart[data_long$trans<3] <- data_long$Tstart[data_long$trans<3] +data_long$age[data_long$trans<3]
data_long$time <- data_long$Tstop-data_long$Tstart

meanlog <- mean(log(data_raw$cov1))
sdlog <- sd(log(data_raw$cov1))


############################
# Data Simulation with hesim
############################

n_pats <- 1000
model_type <- "forward"
model <- load("ground_truthMM.RData")

n_sim <- 100
seeds <- c(555:654)

followup <- 20
dataset <- vector("list", length(n_sim))

start_time <- Sys.time()
dataset <- mclapply(1:n_sim, function(i){
  run_simulation(n_pats, fits_wei, model_type, seeds[i], followup, meanlog, sdlog)
}, mc.cores= detectCores())
end_time <- Sys.time()
elapsed_time <- end_time - start_time
elapsed_time


#setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/Simulation/saved_data")
for (i in 1:n_sim) {
  filename <- paste0("simulation100K_MM_seed_", i, ".RData")
  data <- dataset[[i]]
  save(data, file = filename)
}

