wrapper_functions_SM <- function(data,n_pats,seed,cores_imp){
  
  comp_time <- numeric()
  setwd(here())

  if (n_pats==500){
    if (scheme==2){
      model_dir <- paste0("wrapper_SM/results_500/saved_models_scheme2/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==3){
      model_dir <- paste0("wrapper_SM/results_500/saved_models_scheme3/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==4){
      model_dir <- paste0("wrapper_SM/results_500/saved_models_scheme4/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  
    } else if (scheme==5){
      model_dir <- paste0("wrapper_SM/results_500/saved_models_scheme5/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  }
  } else if (n_pats==2000){
    if (scheme==2){
      model_dir <- paste0("wrapper_SM/results_2K/saved_models_scheme2/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==3){
      model_dir <- paste0("wrapper_SM/results_2K/saved_models_scheme3/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==4){
      model_dir <- paste0("wrapper_SM/results_2K/saved_models_scheme4/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  
    } else if (scheme==5){
      model_dir <- paste0("wrapper_SM/results_2K/saved_models_scheme5/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  }
  } else if (n_pats==5000){
    if (scheme==2){
      model_dir <- paste0("wrapper_SM/results_5K/saved_models_scheme2/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==3){
      model_dir <- paste0("wrapper_SM/results_5K/saved_models_scheme3/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==4){
      model_dir <- paste0("wrapper_SM/results_5K/saved_models_scheme4/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  
    } else if (scheme==5){
      model_dir <- paste0("wrapper_SM/results_5K/saved_models_scheme5/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  }
  } else if (n_pats==10000){
    if (scheme==2){
      model_dir <- paste0("wrapper_SM/results_10K/saved_models_scheme2/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==3){
      model_dir <- paste0("wrapper_SM/results_10K/saved_models_scheme3/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==4){
      model_dir <- paste0("wrapper_SM/results_10K/saved_models_scheme4/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  
    } else if (scheme==5){
      model_dir <- paste0("wrapper_SM/results_10K/saved_models_scheme5/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  }
  }

  n_pats <- length(unique(data$patient_id))
  
  #####################
  # coxph model
  #####################
  
  temp <-  prepare_coxph(data, n_pats)
  
  error <- F
  
  time_coxph <- system.time({
    tryCatch({
      model_cox<- coxph(Surv(Tstart, Tstop ,status) ~ cov1.1 + cov2.1 + cov3.1 + cov1.2 + cov2.2 + cov3.2 + cov1.3 + cov2.3 + cov3.3 + strata(trans), data = temp, method="breslow")
    },
    error = function(e) {
      print(paste("Error during model fitting:", e$message))
      error <<- TRUE
    })
  })[3]
  
  if (error) {
    print(paste("No cox convergence for seed:", seed))
    model_cox <- NULL
  } else {
    print("Model fitted successfully.")
  }
  
  
  comp_time[1] <- as.numeric(round(time_coxph,3))
  
  if (!is.null(model_cox)) {
    save(model_cox, file = file.path(model_dir, "cox_model.RData"))
  } else {
    print("Model Cox is NULL; not saving.")
  }
  
  gc()
  
 
  
  ######################
  # flexsurv model
  ######################
  
  temp <- prepare_flex(data, n_pats)
  fits_gompertz <- vector(mode = "list", length = 3)

  error <- F
  
  time_gomp <- system.time({
    tryCatch({
      for (i in 1:2){
        fits_gompertz[[i]] <- flexsurvreg(Surv(Tstart, Tstop, status) ~ cov1 + cov2 + cov3,
                                        data = subset(temp, trans == i),
                                        dist = "gompertz")
        }
      fits_gompertz[[3]] <- flexsurvreg(Surv(time, status) ~ cov1 + cov2 + cov3,  
                                        data = subset(temp, trans == 3),
                                        dist = "gompertz")
    },
    error = function(e) {
      print(paste("Error during model fitting:", e$message))
      error <<- TRUE
    })
  })[3]
  
  if (error) {
    print(paste("No gomp convergence for seed:", seed))
    fits_gompertz <- NULL
  } else {
    print("Model fitted successfully.")
  }
  
  
  comp_time[2] <- as.numeric(round(time_gomp,3))
  
  if (!is.null(fits_gompertz)) {
    save(fits_gompertz, file = file.path(model_dir, "flexsurv_model.RData"))
  } else {
    print("Model gomp is NULL; not saving.")
  }
  
  gc()
  
  
  
  ####################
  # imputation model
  ####################
  
  temp <- prepare_imputation(data, n_pats)

  
  m <- 30
  type <- "mix"
  
  error <- F
  
  time_imp <- system.time({
    tryCatch({
      results_imp<- run_imputation(temp[[1]], temp[[2]], m, type, cores_imp)
      avg_parameters <- results_imp[[1]]
      all_fits <- results_imp[[2]]
      
    },
    error = function(e) {
      print(paste("Error during model fitting:", e$message))
      error <<- TRUE
    })
  })[3]
  
  if (error) {
    print(paste("No imp convergence for seed:", seed))
    results_imp <- NULL
  } else {
    print("Model fitted successfully.")
  }
  
  
  comp_time[3] <- as.numeric(round(time_imp,3))
  
  if (!is.null(results_imp)) {
    save(results_imp, file = file.path(model_dir, "results_imp.RData"))
  } else {
    print("Model imp is NULL; not saving.")
  }
  
  gc()

  

  save(comp_time, file = file.path(model_dir, "computational_time.RData"))
  
  cat("models completed for seed:", seed, "\n")
  
  }
