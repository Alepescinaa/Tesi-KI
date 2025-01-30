computing_life_expectancy <- function(n_pats, scheme, seed, convergence, t_start, covs){
  
 
  
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

  tmat <- mstate::transMat(x = list(c(2, 3),c(3),c()), names = c("Dementia-free","Dementia", "Death")) 
  # gt_tls <-  totlos.fs.mine(fits_wei, t_start= t_start,  trans=tmat, newdata = covs, t=120)[1,][1:2]
  
  meanlog <- mean(log(fits_wei[[1]]$data$mml$rate[,2]))
  sdlog <- sd(log(fits_wei[[1]]$data$mml$rate[,2]))
  sim_data <- simulation(100000, fits_wei, meanlog, sdlog, covs)
  sim_data_dis <- sim_data$sim_disease(max_t=120, max_age=120)
  sim_data_dis$transition <- 0 #enter in the study
  sim_data_dis$transition[sim_data_dis$from==2 & sim_data_dis$to==3] <- 1 # out from healthy state (dem)
  sim_data_dis$transition[sim_data_dis$from==2 & sim_data_dis$to==4] <- 1 # out from healthy state (death)
  sim_data_dis$transition[sim_data_dis$from==3 & sim_data_dis$to==4] <- 2 # out from dementia
  sim_data_dis$time <- sim_data_dis$time_stop-sim_data_dis$time_start

  mean_time <- sim_data_dis %>%
    group_by(transition) %>%
    summarise(across(time, mean, na.rm = TRUE))

  gt_tls <- as.numeric(unlist(mean_time[2:3,2]))


 
  # ===============
  # EO dataset
  # ===============
  
  #flexsurv_tls_EO <- totlos.fs.mine(fits_gompertz_EO, t_start=t_start,  trans=tmat, newdata = covs, t=105)[1,][1:2]
  #bias_flexsurv_tls_EO <- (flexsurv_tls_EO-gt_tls)/gt_tls
  
  sim_data <- simulation(100000, fits_gompertz_EO, meanlog, sdlog, covs)
  sim_data_dis <- sim_data$sim_disease(max_t=120, max_age=120)
  sim_data_dis$transition <- 0 #enter in the study
  sim_data_dis$transition[sim_data_dis$from==2 & sim_data_dis$to==3] <- 1 # out from healthy state (dem)
  sim_data_dis$transition[sim_data_dis$from==2 & sim_data_dis$to==4] <- 1 # out from healthy state (death)
  sim_data_dis$transition[sim_data_dis$from==3 & sim_data_dis$to==4] <- 2 # out from dementia
  sim_data_dis$time <- sim_data_dis$time_stop-sim_data_dis$time_start
  
  mean_time <- sim_data_dis %>%
    group_by(transition) %>%
    summarise(across(time, mean, na.rm = TRUE))
  
  flexsurv_tls_EO <- as.numeric(unlist(mean_time[2:3,2]))
  bias_flexsurv_tls_EO<- (flexsurv_tls_EO-gt_tls)/gt_tls
  
  # ============
  # coxph
  # ============
  
  if (convergence$coxph[seed]==2) {
    covs_cox <- as.numeric(covs)
    newdata_cox <- data.frame(
      trans = 1:3,
      cov1.1 = c(covs_cox[1], 0, 0),
      cov2.1 = c(covs_cox[2], 0, 0),
      cov3.1 = c(covs_cox[3], 0, 0),
      cov1.2 = c(0, covs_cox[1], 0),
      cov2.2 = c(0, covs_cox[2], 0),
      cov3.2 = c(0, covs_cox[3], 0),
      cov1.3 = c(0, 0, covs_cox[1]),
      cov2.3 = c(0, 0, covs_cox[2]),
      cov3.3 = c(0, 0, covs_cox[3]),
      strata = 1:3
    )
  
    
    msfit_obj <- msfit(model_cox, newdata = newdata_cox, variance=T, trans=tmat)
    
    cox_trans_prob <- probtrans(msfit_obj, predt=min(msfit_obj$Haz$time))[[1]]
    check_neg <- apply(cox_trans_prob, 1, function(row) any(row < 0))
    cox_trans_prob <- cox_trans_prob[!check_neg,]
    time <- cox_trans_prob[,1]
    cox_trans_prob <- cox_trans_prob[,2:4]
    
    coxph_tls<- numeric(ncol(cox_trans_prob))
    diff_time <- diff(time)
    cox_trans_prob <- cox_trans_prob[1:length(diff_time),]
    
    for (i in 1:ncol(cox_trans_prob)) {
      coxph_tls[i] <- sum(cox_trans_prob[, i] * diff(time))
    }
    coxph_tls <- coxph_tls[1:2]
    bias_coxph_tls <- (coxph_tls-gt_tls)/gt_tls
  } else {
    coxph_tls <- rep(NA,2)
    bias_coxph_tls <- rep(NA,2)
  } 
    
 

  # ==============
  # flexsurv
  # ==============
  
  if (convergence$flexsurv[seed]==2){
  #flexsurv_tls <-  totlos.fs.mine(fits_gompertz, t_start=t_start,  trans=tmat, newdata = covs, t=105)[1,][1:2]
  #bias_flexsurv_tls <- (flexsurv_tls-gt_tls)/gt_tls
  
  sim_data <- simulation(100000, fits_gompertz, meanlog, sdlog, covs)
  sim_data_dis <- sim_data$sim_disease(max_t=120, max_age=120)
  sim_data_dis$transition <- 0 #enter in the study
  sim_data_dis$transition[sim_data_dis$from==2 & sim_data_dis$to==3] <- 1 # out from healthy state (dem)
  sim_data_dis$transition[sim_data_dis$from==2 & sim_data_dis$to==4] <- 1 # out from healthy state (death)
  sim_data_dis$transition[sim_data_dis$from==3 & sim_data_dis$to==4] <- 2 # out from dementia
  sim_data_dis$time <- sim_data_dis$time_stop-sim_data_dis$time_start
  
  mean_time <- sim_data_dis %>%
    group_by(transition) %>%
    summarise(across(time, mean, na.rm = TRUE))
  
  flexsurv_tls <- as.numeric(unlist(mean_time[2:3,2]))
  bias_flexsurv_tls <- (flexsurv_tls-gt_tls)/gt_tls
  } else {
    flexsurv_tls <- rep(NA,2)
    bias_flexsurv_tls <- rep(NA,2)
  }
  
  # ========
  # msm
  # ========
  
  # if(convergence$msm[seed]==2){
  #   time <- seq(t_start,105,by=0.1)
  #   time <- time-t_start
  #   msm_probabilities <- matrix(ncol = 3, nrow = 0)
  #   
  #  
  #   for (i in 1:length(time)) {
  #     temp <- pmatrix.msm(model.msm, t = time[i])[1,]
  #     msm_probabilities <- rbind(msm_probabilities, temp)  
  #   }
  #   
  #   msm_tls<- numeric(ncol(msm_probabilities))
  #   
  #   diff_time <- diff(time)
  #   msm_probabilities <- msm_probabilities[1:length(diff_time),]
  #   
  #   for (i in 1:ncol(msm_probabilities)) {
  #     msm_tls[i] <- sum(msm_probabilities[, i] * diff_time)
  #   }
  #   msm_tls <- msm_tls[1:2]
  # }else{
  #   msm_tls <- rep(NA,2)
  # }
  
  
  # if (convergence$msm[seed]==2){
  #   msm_tls <- totlos.msm(model.msm, fromt=0, tot=105-t_start, covariates = list(covs[1],covs[2],covs[3])) [1:2]
  #   bias_msm_tls <- (msm_tls-gt_tls)/gt_tls
  #   
  #   
  # } else{
  #   msm_tls <- rep(NA,2)
  #   bias_msm_tls <- rep(NA,2)
  # }

  if (convergence$msm[seed]==2){
  sim_data <- simulation_new(100000, model.msm, meanlog, sdlog, covs, ind=2)
  sim_data_dis <- sim_data$sim_disease(max_t=120, max_age=120)
  sim_data_dis$transition <- 0 #enter in the study
  sim_data_dis$transition[sim_data_dis$from==2 & sim_data_dis$to==3] <- 1 # out from healthy state (dem)
  sim_data_dis$transition[sim_data_dis$from==2 & sim_data_dis$to==4] <- 1 # out from healthy state (death)
  sim_data_dis$transition[sim_data_dis$from==3 & sim_data_dis$to==4] <- 2 # out from dementia
  sim_data_dis$time <- sim_data_dis$time_stop-sim_data_dis$time_start
  
  mean_time <- sim_data_dis %>%
    group_by(transition) %>%
    summarise(across(time, mean, na.rm = TRUE))
  
  msm_tls <- as.numeric(unlist(mean_time[2:3,2]))
  bias_msm_tls <- (msm_tls-gt_tls)/gt_tls
} else {
  msm_tls <- rep(NA,2)
  bias_msm_tls <- rep(NA,2)
}


  # ========
  # msm_age
  # ========
  
  # elect doens't work as expected
  # to prepare elect
  # 
  # temp$age <- temp$age-t_start
  # temp$onset_age <- temp$onset_age-t_start
  # temp$death_time <- temp$death_time-t_start
  # temp <- temp %>%
  #   group_by(patient_id) %>%
  #   mutate(bsline = ifelse(row_number() == 1, 1, 0)) %>%
  #   ungroup()
  # baseline_data <- temp[temp$bsline==1,]
  # baseline_data$state <- 1
  
  # if(convergence$msm_age[seed]==2){
  # 
  #   mean_age <- mean(model.msm_age$data$mm.cov[,5])
  #   elect_model <- elect( x = model.msm_age, b.covariates = list( age = mean_age, cov1 = cov_means[2], cov2 = cov_means[3], cov3 = cov_means[4]),
  #                         statedistdata = baseline_data, h = 0.1, age.max = 105-t_start) 
  # 
  #   msm_age_tls <- round(elect_model$pnt,3)
  #   msm_age_tls <- msm_age_tls[1:2]
  #   bias_msm_age_tls <- msm_age_tls-gt_tls
  # 
  # 
  # }else{
  #   msm_age_tls <- rep(NA,2)
  #   bias_msm_age_tls <- rep(NA,2)
  # }
  
  # if(convergence$msm_age[seed]==2){
  #   mean_age <- mean(model.msm_age$data$mm.cov[,5])
  #   msm_age_probabilities_cate <- p.matrix.age(model.msm_age, covariates =  list(cov1 = covs[1], cov2 = covs[2], cov3 = covs[3]), age1 = t_start,
  #                                              age2 = 105, int=0.1, ci="none", cores=4)
  #   msm_age_tls <- get_time(msm_age_probabilities_cate, ci= "none")
  #   msm_age_tls <- msm_age_tls[1:2,2]
  #   msm_age_tls <- as.numeric(unlist(msm_age_tls))
  #   bias_msm_age_tls <- msm_age_tls-gt_tls
  # }else{
  #   msm_age_tls <- rep(NA,2)
  #   bias_msm_age_tls <- rep(NA,2)
  # }

 # treating it just as msm understimates
  
  # if(convergence$msm_age[seed]==2){
  #   time <- seq(t_start,105,by=0.1)
  #   time <- time-t_start
  #   msm_age_probabilities <- matrix(ncol = 3, nrow = 0)
  #   
  #   
  #   for (i in 1:length(time)) {
  #     temp <- pmatrix.msm(model.msm_age, t = time[i])[1,]
  #     msm_age_probabilities <- rbind(msm_age_probabilities, temp)
  #   }
  #   
  #   msm_age_tls<- numeric(ncol(msm_age_probabilities))
  #   
  #   diff_time <- diff(time)
  #   msm_age_probabilities <- msm_age_probabilities[1:length(diff_time),]
  #   
  #   for (i in 1:ncol(msm_age_probabilities)) {
  #     msm_age_tls[i] <- sum(msm_age_probabilities[, i] * diff_time)
  #   }
  #   msm_age_tls <- msm_age_tls[1:2]
  #   bias_msm_age_tls <- (msm_age_tls-gt_tls)/gt_tls
  # }else{
  #   msm_age_tls <- rep(NA,2)
  #   bias_msm_age_tls <- rep(NA,2)
  # }
  
  if(convergence$msm_age[seed]==2){
  sim_data <- simulation_new(100000, model.msm_age, meanlog, sdlog, covs, ind=3)
  sim_data_dis <- sim_data$sim_disease(max_t=120, max_age=120)
  sim_data_dis$transition <- 0 #enter in the study
  sim_data_dis$transition[sim_data_dis$from==2 & sim_data_dis$to==3] <- 1 # out from healthy state (dem)
  sim_data_dis$transition[sim_data_dis$from==2 & sim_data_dis$to==4] <- 1 # out from healthy state (death)
  sim_data_dis$transition[sim_data_dis$from==3 & sim_data_dis$to==4] <- 2 # out from dementia
  sim_data_dis$time <- sim_data_dis$time_stop-sim_data_dis$time_start
  
  mean_time <- sim_data_dis %>%
    group_by(transition) %>%
    summarise(across(time, mean, na.rm = TRUE))
  
  msm_age_tls <- as.numeric(unlist(mean_time[2:3,2]))
  bias_msm_age_tls <- (msm_age_tls-gt_tls)/gt_tls
} else {
  msm_age_tls <- rep(NA,2)
  bias_msm_age_tls <- rep(NA,2)
}
  # ======
  # nhm
  # ======
  
if(convergence$nhm[seed]==2){
  sim_data <- simulation_new(100000, model_nhm, meanlog, sdlog, covs, ind=1)
  sim_data_dis <- sim_data$sim_disease(max_t=120, max_age=120)
  sim_data_dis$transition <- 0 #enter in the study
  sim_data_dis$transition[sim_data_dis$from==2 & sim_data_dis$to==3] <- 1 # out from healthy state (dem)
  sim_data_dis$transition[sim_data_dis$from==2 & sim_data_dis$to==4] <- 1 # out from healthy state (death)
  sim_data_dis$transition[sim_data_dis$from==3 & sim_data_dis$to==4] <- 2 # out from dementia
  sim_data_dis$time <- sim_data_dis$time_stop-sim_data_dis$time_start
  
  mean_time <- sim_data_dis %>%
    group_by(transition) %>%
    summarise(across(time, mean, na.rm = TRUE))
  
  nhm_tls <- as.numeric(unlist(mean_time[2:3,2]))
  bias_nhm_tls <- (nhm_tls-gt_tls)/gt_tls
} else {
  nhm_tls <- rep(NA,2)
  bias_nhm_tls <- rep(NA,2)
}
  
  # if(convergence$nhm[seed]==2){
  #   tcrit <- model_nhm$tcrit
  #   time <- seq(t_start,tcrit-1,by=0.1)
  #   nhm_probabilities <- predict(model_nhm, time0= t_start, times= time, covvalue = c(as.numeric(covs[1]),as.numeric(covs[2]),as.numeric(covs[3])))$probabilities
  #   nhm_tls<- numeric(ncol(nhm_probabilities))
  #   
  #   diff_time <- diff(time)
  #   nhm_probabilities <- nhm_probabilities[1:length(diff_time),]
  #   
  #   for (i in 1:ncol(nhm_probabilities)) {
  #     nhm_tls[i] <- sum(nhm_probabilities[, i] * diff_time)
  #   }
  #   nhm_tls <- nhm_tls[1:2]
  #   bias_nhm_tls <- (nhm_tls-gt_tls)/gt_tls
  # }else{
  #   nhm_tls <- rep(NA,2)
  #   bias_nhm_tls <- rep(NA,2)
  # }
  # 
  


  # ==========
  # imputation
  # ==========
  
  # if(convergence$imputation[seed]==2){
  #   imputation_tls <- matrix(0,3,3)
  #   models_imp <- results_imp[[2]]
  #   m <- length(models_imp)
  #   
  #   for (i in 1:m){
  #     temp <-  totlos.fs.mine(models_imp[[i]], t_start=t_start,  trans=tmat, newdata = covs, t=105)
  #     imputation_tls <- imputation_tls + temp
  #   }
  #   imputation_tls <- imputation_tls/m
  #   imputation_tls <- imputation_tls[1,][1:2]
  #   bias_imputation_tls <- (imputation_tls-gt_tls)/gt_tls
  # }else{
  #   imputation_tls <- rep(NA,2)
  #   bias_imputation_tls <- rep(NA,2)
  # }
  
  if(convergence$imputation[seed]==2){
    imputation_tls <- c(0,0)
    models_imp <- results_imp[[2]]
    m <- length(models_imp)
    
    for (i in 1:m){
      sim_data <- simulation(100000, models_imp[[i]], meanlog, sdlog, covs)
      sim_data_dis <- sim_data$sim_disease(max_t=120, max_age=120)
      sim_data_dis$transition <- 0 #enter in the study
      sim_data_dis$transition[sim_data_dis$from==2 & sim_data_dis$to==3] <- 1 # out from healthy state (dem)
      sim_data_dis$transition[sim_data_dis$from==2 & sim_data_dis$to==4] <- 1 # out from healthy state (death)
      sim_data_dis$transition[sim_data_dis$from==3 & sim_data_dis$to==4] <- 2 # out from dementia
      sim_data_dis$time <- sim_data_dis$time_stop-sim_data_dis$time_start
      
      mean_time <- sim_data_dis %>%
        group_by(transition) %>%
        summarise(across(time, mean, na.rm = TRUE))
      
      temp <- as.numeric(unlist(mean_time[2:3,2]))  
      imputation_tls <- imputation_tls + temp
    }
    imputation_tls <- imputation_tls/m
    bias_imputation_tls <- (imputation_tls-gt_tls)/gt_tls
  }else{
    imputation_tls <- rep(NA,2)
    bias_imputation_tls <- rep(NA,2)
  }
  
  # ============
  # wrap up
  # ============
  
  # check quantities
  
  lfe_estimates <- rbind(
    cbind(lfe = flexsurv_tls_EO, model = "flexsurv_EO", seed=seed),
    cbind(lfe = coxph_tls, model = "coxph", seed = seed),
    cbind(lfe = flexsurv_tls, model = "flexsurv", seed = seed),
    cbind(lfe = msm_tls, model = "msm", seed = seed),
    cbind(lfe = msm_age_tls, model = "msm_age", seed = seed),
    cbind(lfe = nhm_tls, model = "nhm", seed = seed),
    cbind(lfe = imputation_tls, model = "imputation", seed = seed) 
  )
  
  lfe_bias <- rbind(
    cbind(lfe = bias_flexsurv_tls_EO, model = "flexsurv_EO", seed=seed),
    cbind(lfe = bias_coxph_tls, model = "coxph", seed = seed),
    cbind(lfe = bias_flexsurv_tls, model = "flexsurv", seed = seed),
    cbind(lfe = bias_msm_tls, model = "msm", seed = seed),
    cbind(lfe = bias_msm_age_tls, model = "msm_age", seed = seed),
    cbind(lfe = bias_nhm_tls, model = "nhm", seed = seed),
    cbind(lfe = bias_imputation_tls, model = "imputation", seed = seed) 
  )
  
  print(seed)
  return(list(lfe_estimates, lfe_bias))
  
}
 
