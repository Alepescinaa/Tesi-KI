run_baseline_bias <- function(n_pats, scheme, seed, convergence){
  
  main_directory <- here()
  
  if (n_pats == 500){
    scheme_dir <- file.path(main_directory, "wrapper_MM/results_500/saved_models_scheme1")
    seed_dir <- file.path(scheme_dir, paste0("seed_", seed))
  } else if (n_pats == 2000){
    scheme_dir <- file.path(main_directory, "wrapper_MM/results_2K/saved_models_scheme1")
    seed_dir <- file.path(scheme_dir, paste0("seed_", seed))
  }else if (n_pats == 5000){
    scheme_dir <- file.path(main_directory, "wrapper_MM/results_5K/saved_models_scheme1")
    seed_dir <- file.path(scheme_dir, paste0("seed_", seed))
  }else if (n_pats == 10000){
    scheme_dir <- file.path(main_directory, "wrapper_MM/results_10K/saved_models_scheme1")
    seed_dir <- file.path(scheme_dir, paste0("seed_", seed))
  }
  
  setwd(seed_dir)
  file <- "fits_gompertz_EO.RData"
  if (file.exists(file)) {
    load(file)
  } else {
    warning(paste("File does not exist:", file))
  }
  
  main_directory <- here()
  
  if (n_pats == 500){
    if (scheme == 2){
      scheme_dir <- file.path(main_directory, "wrapper_MM/results_500/saved_models_scheme2")
    } else if (scheme == 3){
      scheme_dir <- file.path(main_directory, "wrapper_MM/results_500/saved_models_scheme3")
    } else if (scheme == 4){
      scheme_dir <- file.path(main_directory, "wrapper_MM/results_500/saved_models_scheme4")
    } else if (scheme == 5){
      scheme_dir <- file.path(main_directory, "wrapper_MM/results_500/saved_models_scheme5")
    }
  } else if (n_pats == 2000){
    if (scheme == 2){
      scheme_dir <- file.path(main_directory, "wrapper_MM/results_2K/saved_models_scheme2")
    } else if (scheme == 3){
      scheme_dir <-  file.path(main_directory, "wrapper_MM/results_2K/saved_models_scheme3")
    } else if (scheme == 4){
      scheme_dir <-  file.path(main_directory, "wrapper_MM/results_2K/saved_models_scheme4")
    } else if (scheme == 5){
      scheme_dir <-  file.path(main_directory, "wrapper_MM/results_2K/saved_models_scheme5")
    }
  } else if (n_pats == 5000){
    if (scheme == 2){
      scheme_dir <- file.path(main_directory, "wrapper_MM/results_5K/saved_models_scheme2")
    } else if (scheme == 3){
      scheme_dir <- file.path(main_directory, "wrapper_MM/results_5K/saved_models_scheme3")
    } else if (scheme == 4){
      scheme_dir <- file.path(main_directory, "wrapper_MM/results_5K/saved_models_scheme4")
    } else if (scheme == 5){
      scheme_dir <- file.path(main_directory, "wrapper_MM/results_5K/saved_models_scheme5")
    }
  } else if (n_pats == 10000){
    if (scheme == 2){
      scheme_dir <- file.path(main_directory, "wrapper_MM/results_10K/saved_models_scheme2")
    } else if (scheme == 3){
      scheme_dir <- file.path(main_directory, "wrapper_MM/results_10K/saved_models_scheme3")
    } else if (scheme == 4){
      scheme_dir <- file.path(main_directory, "wrapper_MM/results_10K/saved_models_scheme4")
    } else if (scheme == 5){
      scheme_dir <- file.path(main_directory, "wrapper_MM/results_10K/saved_models_scheme5")
    }
  }
  seed_dir <- file.path(scheme_dir, paste0("seed_", seed))
  
  if (dir.exists(seed_dir)) {
    setwd(seed_dir)
    
    files_to_load <- c( 
                       "flexsurv_model.RData", 
                       "model_msm_age.RData", 
                       "model_nhm.RData", 
                       "results_imp.RData")
    for (file in files_to_load) {
      if (file.exists(file)) {
        load(file)
      } else {
        print(file)
        print(seed)
        warning(paste("File does not exist:", file))
        model_nhm <- NULL
      }
    }
  } else {
    warning(paste("Seed directory does not exist:", seed_dir))
  }
  
  # ===========
  # EO dataset
  # ===========
  
  params_EO <- matrix(0, nrow = 3, ncol = 2)
  param_names <- names(fits_gompertz_EO[[1]]$coefficients[1:2])
  colnames(params_EO) <- param_names
  
  for (i in 1:3) {
    for (j in 1:2) {
      params_EO[i, j] <- fits_gompertz_EO[[i]]$coefficients[j]
    }
  }
  
  params_EO <- params_EO[, c(2, 1)]
  
  bias_EO <- (compute_bias(params_EO, ground_truth_params[,1:2]))/ground_truth_params[,1:2]
  
  
  
   # ============
   # flexsurv
   # ============
  
  if (convergence$flexsurv[seed] == 2) {
    params_flexsurv <- matrix(0, nrow = 3, ncol = 2)
    param_names <- names(fits_gompertz[[1]]$coefficients[1:2])
    colnames(params_flexsurv) <- param_names
    
    for (i in 1:3) {
      for (j in 1:2) {
        params_flexsurv[i, j] <- fits_gompertz[[i]]$coefficients[j]
      }
    }
    
    params_flexsurv <- params_flexsurv[, c(2, 1)]
  
    bias_flexsurv <- compute_bias(params_flexsurv, ground_truth_params[,1:2])/ground_truth_params[,1:2]
  } else {
    params_flexsurv <- matrix(NA, nrow = 3, ncol = 2)
    bias_flexsurv <- matrix(NA, nrow = 3, ncol = 2)
    colnames(bias_flexsurv) <- colnames(ground_truth_params[,1:2])
    rownames(bias_flexsurv) <- rownames(ground_truth_params[,1:2])
  }
  gc()
  

   # ============
   # msm + age
   # ============
  
  if (convergence$msm_age[seed]==2){
  
    min_age <- min(model.msm_age$data[[1]]$age)
    max_age <- max(model.msm_age$data[[1]]$age)
    #mean_age <- mean(model.msm_age$data[[1]]$age)
    covs <- colMeans(model.msm_age$data[[1]][,1:3])
    haz <- hazards_mine(model.msm_age, b.covariates = list(age = min_age , cov1 = covs[1], cov2 = covs[2], cov3 = covs[3]), no.years = max_age-min_age, CI = F)
    
    # Assuming this hazards come from fitting a gompertz model I wanna retrieve for each transition
    # shape and rate value, parameters of the distribution
    # h(t)=rate*exp(shape*t)
    # log(h(t))= log(rate) + shape*t
    # can be seen as y(t)= a+b*t
    shape <- numeric()
    rate <- numeric()
    
    for (i in 1:3){
      age_grid <- haz[[2]]
      y <- log(as.numeric(unlist(haz[[1]][[i]])))
      reg_model <- lm(y ~ age_grid)
      rate[[i]] <- reg_model$coefficients[1]
      shape[[i]]<- reg_model$coefficients[2]
    }
    
    params_msm_age <- cbind(rate,shape)
    bias_msm_age <- compute_bias(params_msm_age, ground_truth_params[,1:2])/ground_truth_params[,1:2]
    
  } else {
    params_msm_age <- matrix(NA, nrow = 3, ncol = 2)
    bias_msm_age <- matrix(NA, nrow = 3, ncol = 2)
    colnames(bias_msm_age) <- colnames(ground_truth_params[,1:2])
    rownames(bias_msm_age) <- rownames(ground_truth_params[,1:2])
  }
  
  gc()
  
  # ============
  # nhm
  # ============
  
  if (convergence$nhm[seed]==2){
    params_nhm <- matrix(model_nhm$par[1:6], nrow = 3, ncol = 2)
    colnames(params_nhm) <- colnames(ground_truth_params[,1:2])
    
    bias_nhm <- compute_bias(params_nhm, ground_truth_params[,1:2])/ground_truth_params[,1:2]
  } else{
    params_nhm <- matrix(NA, nrow = 3, ncol = 2)
    bias_nhm <- matrix(NA, nrow = 3, ncol = 2)
    colnames(bias_nhm) <- colnames(ground_truth_params[,1:2])
    rownames(bias_nhm) <- rownames(ground_truth_params[,1:2])
    
  }
  
  gc()
   
  
   # ============
   # imputation
   # ============

  if (convergence$imputation[seed]==2){
    params_imp <- results_imp[[1]]
    params_imp <- params_imp[, c(2, 1)]
    
    bias_imputation <- compute_bias(params_imp, ground_truth_params[,1:2])/ground_truth_params[,1:2]
  } else{
    bias_imputation <- matrix(NA, nrow = 3, ncol = 2)
    colnames(bias_imputation) <- colnames(ground_truth_params[,1:2])
    rownames(bias_imputation) <- rownames(ground_truth_params[,1:2])
    
  }
  
  
  bias_tot <- rbind(
    cbind(bias_EO, model = "flexsurv_EO", seed=seed, transition = c(1,2,3)),
    cbind(bias_flexsurv, model = "flexsurv", seed = seed, transition = c(1,2,3)),
    cbind(bias_msm_age, model = "msm_age", seed = seed, transition = c(1,2,3)),
    cbind(bias_nhm, model = "nhm", seed = seed, transition = c(1,2,3)),
    cbind(bias_imputation, model = "imputation", seed = seed, transition = c(1,2,3))
  )
  

  estimates <- rbind(
    cbind(params_EO, model = "flexsurv_EO", seed=seed, transition = c(1,2,3)),
    cbind(params_flexsurv, model = "flexsurv", seed = seed, transition = c(1,2,3)),
    cbind(params_msm_age, model = "msm_age", seed = seed, transition = c(1,2,3)),
    cbind(params_nhm, model = "nhm", seed = seed, transition = c(1,2,3)),
    cbind(params_imp, model = "imputation", seed = seed, transition = c(1,2,3))
  )
  
  
  return (list(bias_tot, estimates))
}

