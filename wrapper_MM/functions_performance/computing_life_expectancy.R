computing_life_expectancy <- function(){
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
    
    files_to_load <- c("cox_model.RData", 
                       "flexsurv_model.RData", 
                       "msm_model.RData", 
                       "model_msm_age.RData", 
                       "model_nhm.RData", 
                       "results_imp.RData", 
                       "computational_time.RData")
    for (file in files_to_load) {
      if (file.exists(file)) {
        load(file)
      } else {
        warning(paste("File does not exist:", file))
        model_nhm <- NULL
      }
    }
  } else {
    warning(paste("Seed directory does not exist:", seed_dir))
  }
  
  # totlos.fs computes expected amount of time spent in state s for a time-inhomogeneous,
  # continuous-time Markov multi-state process that starts in state r up to a maximum time t.
  # This is defined as the integral of the corresponding transition probability up to that time.
  # The argument x is either a model fitted with flexsurv or a list (one for each transition)
  # We are interested in first line of totlos.fs since we know the starting state is 1 (totlos.msm directly computes from state 1)
  
  # ============================
  # life expectancy ground truth
  # ============================
  
  generator_model <- fits_gompertz_EO
  gt <- ground_truth_params[,c(2,1,3,4,5)]
  generator_model[[1]]$coefficients <- gt[1,]
  generator_model[[2]]$coefficients <- gt[2,]
  generator_model[[3]]$coefficients <- gt[3,]
  
  generator_model[[1]]$res.t[,1] <-  gt[1,]
  generator_model[[2]]$res.t[,1] <-  gt[2,]
  generator_model[[3]]$res.t[,1] <-  gt[3,]
  
  generator_model[[1]]$res[,1] <-  gt[1,]
  generator_model[[2]]$res[,1] <-  gt[2,]
  generator_model[[3]]$res[,1] <-  gt[3,]
  
  generator_model[[1]]$res[2,1] <- exp(gt[1,2])
  generator_model[[2]]$res[2,2] <- exp(gt[2,2])
  generator_model[[3]]$res[2,3] <- exp(gt[3,2])
  
  cov_means <- colMeans(generator_model[[1]]$data$mml$rate)
  
  newdata <- data.frame(
    cov1 = cov_means[2], 
    cov2 = cov_means[3],
    cov3 = cov_means[4]
  )
  
  tmat <- mstate::transMat(x = list(c(2, 3),c(3),c()), names = c("Dementia-free","Dementia", "Death")) 
  totlos.fs(generator_model, trans=tmat, newdata = newdata, t=100)
  
  # ===============
  # EO dataset
  # ===============

  totlos.fs(fits_gompertz_EO, trans=tmat, newdata = newdata, t=100)
  
  
  # ============
  # coxph
  # ============
  
  time <- seq(0, 100, by = 0.1)
  
  p_no_dem <- predict(model_cox, type = "expected" ) 
  surv_fit <- survfit(model_cox[[1]], newdata = cov_means, times = time)
  
  cox_tls <- numeric(length(surv_fit$strata))
  
  # Compute total length of stay for each status
  for (i in seq_along(cox_tls)) {
    # Calculate the survival probabilities for the i-th state
    survival_probabilities <- surv_fit$surv[i, ]  # Extract survival probabilities for the i-th state
    cox_tls[i] <- sum(survival_probabilities * diff(c(0, time)))
  }
  
  # Print total length of stay for each status
  names(cox_tls) <- names(surv_fit$strata)  # Assign names based on strata
  print(cox_tls)
  
  
  # ==============
  # flexsurv
  # ==============
  
  totlos.fs(fits_gompertz, trans=tmat, newdata = newdata, t=100)
  
  # ========
  # msm
  # ========
  
  totlos.msm(model.msm, tot=Inf)  #check it by hand

  # ========
  # msm_age
  # ========
  
  totlos.msm(model.msm_age, tot=Inf)  #check elect
  
  # ======
  # nhm
  # ======
  
  time <- seq(0,100,by=0.1)
  nhm_probabilities <- predict(model_nhm, times= time)$probabilities # automatically uses means of covs
  nhm_tls<- numeric(ncol(nhm_probabilities))
  
  for (i in 1:ncol(nhm_probabilities)) {
    nhm_tls[i] <- sum(nhm_probabilities[, i] * diff(c(0, time)))
  }
  
  nhm_tls
  
  # ==========
  # imputation
  # ==========
  
  imputation_tls <- matrix(0,3,3)
  models_imp <- results_imp[[2]]
  m <- length(models_imp)
  
  for (i in 1:m){
    temp <- totlos.fs(models_imp[[i]], trans=tmat, newdata = newdata, t=30)
    imputation_tls <- imputation_tls + temp
  }
  
  imputation_tls <- imputation_tls/m
}

