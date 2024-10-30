gt_flexsurv <- function(n_pats, seed ){
  
  setwd(here())
  load("wrapper_MM/ground_truthMM.RData")
 
  
  if (n_pats == 500){
    model_dir <- paste0("wrapper_MM/results_500/saved_models_scheme1/seed_", seed)
    dir.create(model_dir, showWarnings = FALSE, recursive= T)
    load("Simulated_data_MM/simulation500_MM_all.RData")
    EO_data <- dataset_all_MM_500[[seed]][[1]]
  } else if (n_pats == 2000){
    model_dir <- paste0("wrapper_MM/results_2K/saved_models_scheme1/seed_", seed)
    dir.create(model_dir, showWarnings = FALSE, recursive= T)
    load("Simulated_data_MM/simulation2K_MM_all.RData")
    EO_data <- dataset_all_MM_2K[[seed]][[1]]
  }else if (n_pats == 5000){
    model_dir <- paste0("wrapper_MM/results_5K/saved_models_scheme1/seed_", seed)
    dir.create(model_dir, showWarnings = FALSE, recursive= T)
    load("Simulated_data_MM/simulation5K_MM_all.RData")
    EO_data <- dataset_all_MM_2K[[seed]][[1]]
  }else if (n_pats == 10000){
    model_dir <- paste0("wrapper_MM/results_10K/saved_models_scheme1/seed_", seed)
    dir.create(model_dir, showWarnings = FALSE, recursive= T)
    load("Simulated_data_MM/simulation10K_MM_all.RData")
    EO_data <- dataset_all_MM_10K[[seed]][[1]]
  }
  
  tmat <- mstate::transMat(x = list(c(2, 3),c(3),c()), names = c("Dementia-free","Dementia", "Death")) 
  
  data_long_EO <- msprep(data = EO_data, trans = tmat, 
                         time = c(NA, "onset_age", "death_time"), 
                         status = c(NA, "onset", "dead"), 
                         keep = c("age","cov1", "cov2", "cov3"),
                         id="patient_id")
  
  data_long_EO$Tstart[data_long_EO$trans<3] <- data_long_EO$Tstart[data_long_EO$trans<3] +data_long_EO$age[data_long_EO$trans<3]
  
  data_long_EO$time <- data_long_EO$Tstop-data_long_EO$Tstart
  
  n_trans <- max(tmat, na.rm = TRUE)
  fits_gompertz_EO <- vector(mode = "list", length = n_trans)
  
  for (i in 1:3) {
    fits_gompertz_EO[[i]] <- flexsurvreg(Surv(Tstart, Tstop, status) ~ cov1 + cov2 + cov3, 
                                    data = subset(data_long_EO, trans == i), 
                                    dist = "gompertz")
  }
  
  save(fits_gompertz_EO, file = file.path(model_dir, "fits_gompertz_EO.RData"))
}

