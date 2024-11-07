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
  
  setwd(here())
  
  if (n_pats == 500){
    load("./Simulated_data_MM/simulation500_MM_all.RData")
    data <- dataset_all_MM_500
    }else if(n_pats == 2000){
    load("./Simulated_data_MM/simulation2K_MM_all.RData")
    data <- dataset_all_MM_2K
  } else if(n_pats == 5000){
    load("./Simulated_data_MM/simulation5K_MM_all.RData")
    data <- dataset_all_MM_5K
  } else if (n_pats == 10000){
    load("./Simulated_data_MM/simulation10K_MM_all.RData")
    data <- dataset_all_MM_10K}
  
  data <- data[[seed]][[scheme]]
  t_start <- min(data$age)

  
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
 
  gt_tls <-  totlos.fs.mine(generator_model, t_start= t_start,  trans=tmat, newdata = newdata, t=105)[1,]
 
  # ===============
  # EO dataset
  # ===============
  
  flexsurv_tls_EO <- totlos.fs.mine(fits_gompertz_EO, t_start=60,  trans=tmat, newdata = newdata, t=105)[1,]
  
  # ============
  # coxph
  # ============
  
  newdata_cox <- data.frame(
    trans = 1:3,
    cov1.1 = c(cov_means[2], 0, 0),
    cov2.1 = c(cov_means[3], 0, 0),
    cov3.1 = c(cov_means[4], 0, 0),
    cov1.2 = c(0, cov_means[2], 0),
    cov2.2 = c(0, cov_means[3], 0),
    cov3.2 = c(0, cov_means[4], 0),
    cov1.3 = c(0, 0, cov_means[2]),
    cov2.3 = c(0, 0, cov_means[3]),
    cov3.3 = c(0, 0, cov_means[4]),
    strata = 1:3
  )
  
  msfit_obj <- msfit(model_cox, newdata = newdata_cox, variance=T, trans=tmat)

  cox_trans_prob <- probtrans(msfit_obj, predt=t_start)[[1]]
  check_neg <- apply(cox_trans_prob, 1, function(row) any(row < 0))
  cox_trans_prob <- cox_trans_prob[!check_neg,]
  time <- cox_trans_prob[,1]
  cox_trans_prob <- cox_trans_prob[,2:4]
  
  coxph_tls<- numeric(ncol(cox_trans_prob))
  
  for (i in 1:ncol(cox_trans_prob)) {
    coxph_tls[i] <- sum(cox_trans_prob[, i] * diff(time))
  }
  
  # ==============
  # flexsurv
  # ==============
  
  flexsurv_tls <-  totlos.fs.mine(fits_gompertz, t_start=60,  trans=tmat, newdata = newdata, t=105)[1,]
  
  # ========
  # msm
  # ========
  
  msm_tls <- totlos.msm(model.msm, fromt=t_start, tot=105)  

  # ========
  # msm_age
  # ========
  

  temp <- data
  temp$age <- temp$age - t_start
  temp <- temp %>%
    group_by(patient_id) %>%
    mutate(bsline = ifelse(row_number() == 1, 1, 0)) %>%
    ungroup()
  
  baseline_data <- temp[temp$bsline==1,]
  baseline_data$state <- 1
  
  max_age <- max(temp$age)
  mean_age <- mean(temp$age) 
  
  elect_model <- elect( x = model.msm_age, b.covariates = list( age = mean_age, cov1 = cov_means[2], cov2 = cov_means[3], cov3 = cov_means[4]),
                        statedistdata = baseline_data, h = 0.1, age.max = max_age)
  
  msm_age_tls <- round(elect_model$pnt,3)
  
  # ======
  # nhm
  # ======
  
  time <- seq(60,105,by=0.01)
  nhm_probabilities <- predict(model_nhm, times= time)$probabilities # automatically uses means of covs
  nhm_tls<- numeric(ncol(nhm_probabilities))
  
  for (i in 1:ncol(nhm_probabilities)) {
    nhm_tls[i] <- sum(nhm_probabilities[, i] * diff(time))
  }
  
  
  # ==========
  # imputation
  # ==========
  
  imputation_tls <- matrix(0,3,3)
  models_imp <- results_imp[[2]]
  m <- length(models_imp)
  
  for (i in 1:m){
    temp <-  totlos.fs.mine(models_imp[[i]], t_start=60,  trans=tmat, newdata = newdata, t=105)
    imputation_tls <- imputation_tls + temp
  }
  imputation_tls <- imputation_tls/m
  imputation_tls <- imputation_tls[1,]
  
}
 
