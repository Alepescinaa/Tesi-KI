library(dplyr)

#######################################################################

# setwd("~/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code")
# 
# dataset_all_MM <-  vector("list", 100)
# for (i in 1:100) {
#   file_name <- paste0("./Simulated_data_MM/simulation100K_MM_seed_", i, ".RData")
#   load(file_name)
#   dataset_all_MM[[i]] <- datai
# }
# save(dataset_all_MM, file= "simulation100K_MM_all.RData")

# dataset_all_SM <-  vector("list", 100)
# for (i in 1:100) {
#   file_name <- paste0("./Simulated_data_SM/simulation100K_SM_seed_", i, ".RData")
#   load(file_name)
#   dataset_all_SM[[i]] <- datai
# }
# save(dataset_all_SM, file= "simulation100K_SM_all.RData")
# 

load("./Simulated_data_MM/simulation100K_MM_all.Rdata")

#######################################################################
dataset_all_MM_500 <- dataset_all_MM
n_pats = 500

for (i in 1:100){
  for (j in 1:5){
    data_tot <- dataset_all_MM_500[[i]]
    data <- data_tot[[1]]
    tiny_data <- data %>%
      group_by(onset,dead) %>%
      sample_n(size = round(n_pats * n()/nrow(data))) %>%
      ungroup()
    keep_id <- tiny_data$patient_id
    temp <- dataset_all_MM_500[[i]][[j]]
    dataset_all_MM_500[[i]][[j]] <- temp[temp$patient_id %in% keep_id, ]
  }}
save(dataset_all_MM_500, file="simulation500_MM_all.Rdata")

#######################################################################
dataset_all_MM_5K <- dataset_all_MM
n_pats = 5000

for (i in 1:100){
  for (j in 1:5){
    data_tot <- dataset_all_MM_5K[[i]]
    data <- data_tot[[1]]
    tiny_data <- data %>%
      group_by(onset,dead) %>%
      sample_n(size = round(n_pats * n()/nrow(data))) %>%
      ungroup()
    keep_id <- tiny_data$patient_id
    temp <- dataset_all_MM_5K[[i]][[j]]
    dataset_all_MM_5K[[i]][[j]] <- temp[temp$patient_id %in% keep_id, ]
    }}
save(dataset_all_MM_5K, file="simulation5K_MM_all.Rdata")

#######################################################################
dataset_all_MM_2K <- dataset_all_MM
n_pats = 2000

for (i in 1:100){
  for (j in 1:5){
    data_tot <- dataset_all_MM_2K[[i]]
    data <- data_tot[[1]]
    tiny_data <- data %>%
      group_by(onset,dead) %>%
      sample_n(size = round(n_pats * n()/nrow(data))) %>%
      ungroup()
    keep_id <- tiny_data$patient_id
    temp <- dataset_all_MM_2K[[i]][[j]]
    dataset_all_MM_2K[[i]][[j]] <- temp[temp$patient_id %in% keep_id, ]
  }}
save(dataset_all_MM_2K, file="simulation2K_MM_all.Rdata")

#######################################################################

dataset_all_MM_10K <- dataset_all_MM
n_pats = 10000

for (i in 1:100){
  for (j in 1:5){
    data_tot <- dataset_all_MM_10K[[i]]
    data <- data_tot[[1]]
    tiny_data <- data %>%
      group_by(onset,dead) %>%
      sample_n(size = round(n_pats * n()/nrow(data))) %>%
      ungroup()
    keep_id <- tiny_data$patient_id
    temp <- dataset_all_MM_10K[[i]][[j]]
    dataset_all_MM_10K[[i]][[j]] <- temp[temp$patient_id %in% keep_id, ]
  }}
save(dataset_all_MM_10K, file="simulation10K_MM_all.Rdata")

#######################################################################

load("./Simulated_data_SM/simulation100K_SM_all.Rdata")

#######################################################################
dataset_all_SM_5K <- dataset_all_SM
n_pats = 5000

for (i in 1:100){
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
save(dataset_all_SM_5K, file="simulation5K_SM_all.Rdata")

#######################################################################
dataset_all_SM_2K <- dataset_all_SM
n_pats = 2000

for (i in 1:100){
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
save(dataset_all_SM_2K, file="simulation2K_SM_all.Rdata")

#######################################################################

dataset_all_SM_10K <- dataset_all_SM
n_pats = 10000

for (i in 1:100){
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
save(dataset_all_SM_10K, file="simulation10K_SM_all.Rdata")
