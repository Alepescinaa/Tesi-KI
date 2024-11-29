
run_performance_bias <- function(n_pats, scheme, seed, convergence){
  
  main_directory <- here()

  if (n_pats == 500){
    scheme_dir <- file.path(main_directory, "wrapper_SM/results_500/saved_models_scheme1")
    seed_dir <- file.path(scheme_dir, paste0("seed_", seed))
  } else if (n_pats == 2000){
    scheme_dir <- file.path(main_directory, "wrapper_SM/results_2K/saved_models_scheme1")
    seed_dir <- file.path(scheme_dir, paste0("seed_", seed))
  }else if (n_pats == 5000){
    scheme_dir <- file.path(main_directory, "wrapper_SM/results_5K/saved_models_scheme1")
    seed_dir <- file.path(scheme_dir, paste0("seed_", seed))
  }else if (n_pats == 10000){
    scheme_dir <- file.path(main_directory, "wrapper_SM/results_10K/saved_models_scheme1")
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
      scheme_dir <- file.path(main_directory, "wrapper_SM/results_500/saved_models_scheme2")
    } else if (scheme == 3){
      scheme_dir <- file.path(main_directory, "wrapper_SM/results_500/saved_models_scheme3")
    } else if (scheme == 4){
      scheme_dir <- file.path(main_directory, "wrapper_SM/results_500/saved_models_scheme4")
    } else if (scheme == 5){
      scheme_dir <- file.path(main_directory, "wrapper_SM/results_500/saved_models_scheme5")
    }
  } else if (n_pats == 2000){
    if (scheme == 2){
      scheme_dir <- file.path(main_directory, "wrapper_SM/results_2K/saved_models_scheme2")
    } else if (scheme == 3){
      scheme_dir <-  file.path(main_directory, "wrapper_SM/results_2K/saved_models_scheme3")
    } else if (scheme == 4){
      scheme_dir <-  file.path(main_directory, "wrapper_SM/results_2K/saved_models_scheme4")
    } else if (scheme == 5){
      scheme_dir <-  file.path(main_directory, "wrapper_SM/results_2K/saved_models_scheme5")
    }
  } else if (n_pats == 5000){
    if (scheme == 2){
      scheme_dir <- file.path(main_directory, "wrapper_SM/results_5K/saved_models_scheme2")
    } else if (scheme == 3){
      scheme_dir <- file.path(main_directory, "wrapper_SM/results_5K/saved_models_scheme3")
    } else if (scheme == 4){
      scheme_dir <- file.path(main_directory, "wrapper_SM/results_5K/saved_models_scheme4")
    } else if (scheme == 5){
      scheme_dir <- file.path(main_directory, "wrapper_SM/results_5K/saved_models_scheme5")
    }
  } else if (n_pats == 10000){
    if (scheme == 2){
      scheme_dir <- file.path(main_directory, "wrapper_SM/results_10K/saved_models_scheme2")
    } else if (scheme == 3){
      scheme_dir <- file.path(main_directory, "wrapper_SM/results_10K/saved_models_scheme3")
    } else if (scheme == 4){
      scheme_dir <- file.path(main_directory, "wrapper_SM/results_10K/saved_models_scheme4")
    } else if (scheme == 5){
      scheme_dir <- file.path(main_directory, "wrapper_SM/results_10K/saved_models_scheme5")
    }
  }
  seed_dir <- file.path(scheme_dir, paste0("seed_", seed))
  
  if (dir.exists(seed_dir)) {
    setwd(seed_dir)
    
    files_to_load <- c("cox_model.RData", 
                       "flexsurv_model.RData", 
                       "model_nhm.RData", 
                       "results_imp.RData", 
                       "computational_time.RData")
    for (file in files_to_load) {
      if (file.exists(file)) {
        load(file)
      } else {
        warning(paste("File does not exist:", file))
        file <- sub("\\.Rdata$", "", file, ignore.case = T) 
        if(file == "cox_model")
          model_cox <- NULL
        if(file == "flexsurv_model")
          fits_gompertz <- NULL
        assign(file,NULL)
        print(seed)
        print(file)
        
      }
    }
    
  } else {
    warning(paste("Seed directory does not exist:", seed_dir))
  }
  
  
  # ===========
  # EO dataset
  # ===========
  
  params_EO <- matrix(0, nrow = 3, ncol = 5)
  param_names <- names(fits_gompertz_EO[[1]]$coefficients)
  colnames(params_EO) <- param_names

  for (i in 1:3){
    for (j in 1:5){
      params_EO[i,j] <- fits_gompertz_EO[[i]] $coefficients[j]
    }
  }

  params_EO <- params_EO[, c(2, 1, 3, 4, 5)]
  params_EO <- cbind(params_EO, exp(params_EO[,3]), exp(params_EO[,4]), exp(params_EO[,4]))
  colnames(params_EO)[6:8] <- c("exp(cov1)", "exp(cov2)", "exp(cov3)")

  bias_EO <- compute_bias(params_EO, ground_truth_params)


 #  # =========
 #  # coxph
 #  # =========
 #  
  if (convergence$coxph[seed]==2) {

    param_names <- c("cov1","cov2", "cov3")
    params_coxph <- matrix(model_cox$coefficients, 3, 3, byrow = T)
    colnames(params_coxph) <- param_names
    params_coxph <- cbind(params_coxph, exp(params_coxph[,1]), exp(params_coxph[,2]), exp(params_coxph[,3]))
    colnames(params_coxph)[4:6] <- c("exp(cov1)", "exp(cov2)", "exp(cov3)")

    bias_coxph <- compute_bias(params_coxph, ground_truth_params)
    bias_coxph[,1:2] <- NA
  }else{
    params_coxph <- matrix(NA, nrow = 3, ncol = ncol(ground_truth_params))
    bias_coxph <- matrix(NA, nrow = 3, ncol = ncol(ground_truth_params))
    colnames(bias_coxph) <- colnames(ground_truth_params)
    rownames(bias_coxph) <- rownames(ground_truth_params)

  }

 #  # ============
 #  # flexsurv
 #  # ============
 #  
  if (convergence$flexsurv[seed]==2){
    params_flexsurv <- matrix(0, nrow = 3, ncol = 5)
    param_names <- names(fits_gompertz[[1]]$coefficients)
    colnames(params_flexsurv) <- param_names

    for (i in 1:3){
      for (j in 1:5){
        params_flexsurv[i,j] <- fits_gompertz[[i]] $coefficients[j]
      }
    }

    params_flexsurv <- params_flexsurv[, c(2, 1, 3, 4, 5)]
    params_flexsurv <- cbind(params_flexsurv, exp(params_flexsurv[,3]), exp(params_flexsurv[,4]), exp(params_flexsurv[,4]))
    colnames(params_flexsurv)[6:8] <- c("exp(cov1)", "exp(cov2)", "exp(cov3)")

    bias_flexsurv <- compute_bias(params_flexsurv, ground_truth_params)

  }else{
    params_flexsurv <- matrix(NA, nrow = 3, ncol = ncol(ground_truth_params))
    bias_flexsurv <- matrix(NA, nrow = 3, ncol = ncol(ground_truth_params))
    colnames(bias_flexsurv) <- colnames(ground_truth_params)
    rownames(bias_flexsurv) <- rownames(ground_truth_params)
  }
  gc()

 #  # ============
 #  # smms
 #  # ============
 #  
 #  if (convergence$msm[seed]==2){
 #    params_msm <- matrix(model.msm$estimates[4:12], nrow = 3, ncol = 3) #1:3 rate
 #    params_msm <- cbind(params_msm, exp(params_msm[,1]), exp(params_msm[,2]), exp(params_msm[,3]))
 #    colnames(params_msm) <- colnames(ground_truth_params)[3:8]
 #    
 #    bias_msm <- compute_bias(params_msm, ground_truth_params)
 #    bias_msm[,1:2] <-NA
 #  } else {
 #    params_msm <- matrix(NA, nrow = 3, ncol = ncol(ground_truth_params))
 #    bias_msm <- matrix(NA, nrow = 3, ncol = ncol(ground_truth_params))
 #    colnames(bias_msm) <- colnames(ground_truth_params)
 #    rownames(bias_msm) <- rownames(ground_truth_params)
 # }
 #  

 #  # ============
 #  # nhm
 #  # ============
 #  
  if (convergence$nhm[seed]==2){
    params_nhm <- matrix(model_nhm$par, nrow = 3, ncol = 5)
    params_nhm <- cbind(params_nhm, exp(params_nhm[,3]), exp(params_nhm[,4]), exp(params_nhm[,5]))
    colnames(params_nhm) <- colnames(ground_truth_params)

    bias_nhm <- compute_bias(params_nhm, ground_truth_params)
  } else{
    params_nhm <- matrix(NA, nrow = 3, ncol = ncol(ground_truth_params))
    bias_nhm <- matrix(NA, nrow = 3, ncol = ncol(ground_truth_params))
    colnames(bias_nhm) <- colnames(ground_truth_params)
    rownames(bias_nhm) <- rownames(ground_truth_params)

  }

  gc()
  #
 #  # ============
 #  # imputation
 #  # ============
 #  
  if (convergence$imputation[seed]==2){
    params_imp <- results_imp[[1]]
    params_imp <- params_imp[, c(2, 1, 3, 4, 5)]
    params_imp <- cbind(params_imp, exp(params_imp[,3]), exp(params_imp[,4]), exp(params_imp[,4]))
    colnames(params_imp)[6:8] <- c("exp(cov1)", "exp(cov2)", "exp(cov3)")

    bias_imputation <- compute_bias(params_imp, ground_truth_params)
  } else{
    bias_imputation <- matrix(NA, nrow = 3, ncol = ncol(ground_truth_params))
    colnames(bias_imputation) <- colnames(ground_truth_params)
    rownames(bias_imputation) <- rownames(ground_truth_params)

   }
 # 

  bias_tot <- rbind(
    cbind(bias_EO, model = "flexsurv_EO", seed=seed, transition = c(1,2,3)),
    cbind(bias_coxph, model = "coxph", seed = seed, transition = c(1,2,3)),
    cbind(bias_flexsurv, model = "flexsurv", seed = seed, transition = c(1,2,3)),
    #cbind(bias_smms, model = "smms", seed = seed, transition = c(1,2,3)),
    cbind(bias_nhm, model = "nhm", seed = seed, transition = c(1,2,3)),
    cbind(bias_imputation, model = "imputation", seed = seed, transition = c(1,2,3))
  )

  rate <- c(NA,NA,NA)
  shape <- c(NA,NA,NA)
  if (convergence$coxph[seed]==2)
    params_coxph <- cbind(rate,shape,params_coxph)
  estimates <- rbind(
    cbind(params_EO, model = "flexsurv_EO", seed=seed, transition = c(1,2,3)),
    cbind(params_coxph, model = "coxph", seed = seed, transition = c(1,2,3)),
    cbind(params_flexsurv, model = "flexsurv", seed = seed, transition = c(1,2,3)),
    #cbind(params_smms, model = "smms", seed = seed, transition = c(1,2,3)),
    cbind(params_nhm, model = "nhm", seed = seed, transition = c(1,2,3)),
    cbind(params_imp, model = "imputation", seed = seed, transition = c(1,2,3))
  )


  return (list(bias_tot, estimates))
}

