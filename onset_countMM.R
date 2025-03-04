library(survival)
library(dplyr)
library(flexsurv)
library(mstate)
library(future.apply)
library(here)


# choose the sample size and upload accordingly the datset, either 500, 2K, 5K
setwd(here())
n_pats <- 2000 # CHANGE HERE


source("./wrapper_MM/functions_wrapper/prepare_imputation.R")
source("./wrapper_MM/functions_wrapper/fit_model.R")
source("./wrapper_MM/functions_wrapper/averaging_params.R")
source("./wrapper_MM/functions_wrapper/run_imputation_dataset.R")
source("./wrapper_MM/functions_wrapper/onset_count.R")

cores <- 4

load("./Simulated_data_MM/simulation2K_MM_all.RData")
results <- list()  # Initialize the results list

for (scheme in 2:5) {  
  results[[scheme]] <- lapply(1:100, function(seed) {  # Use lapply() for sequential execution
    data <- dataset_all_MM_2K[[seed]][[scheme]]  
    n_pats <- length(unique(data$patient_id))  
    gt_onset <- sum(dataset_all_MM_2K[[seed]][[1]]$onset) / n_pats  
    scheme_onset <- data %>%
      group_by(patient_id) %>%
      summarise(dementia_observed = as.integer(any(onset == 1))) %>%
      summarise(mean_onset = mean(dementia_observed))
    imputated_onset <- onset_count(data, n_pats)  
    print(seed)
    return(c(seed, gt_onset, as.numeric(scheme_onset), imputated_onset))  
  })
}

save(results, file = "results_onset.RData")


df_2 <- do.call(rbind, results[[2]])
colnames(df_2) <- c("Seed", "gt_onset", "scheme_onset", "Imp_onset")
colMeans(df_2[,2:4])

df_3 <- do.call(rbind, results[[3]])
colnames(df_3) <- c("Seed", "gt_onset", "scheme_onset", "Imp_onset")
colMeans(df_3[,2:4])

df_4 <- do.call(rbind, results[[4]])
colnames(df_4) <- c("Seed", "gt_onset", "scheme_onset", "Imp_onset")
colMeans(df_4[,2:4])

df_5 <- do.call(rbind, results[[5]])
colnames(df_5) <- c("Seed", "gt_onset", "scheme_onset", "Imp_onset")
colMeans(df_5[,2:4])
