computing_life_expectancy <- function(n_pats, scheme, seed, convergence, t_start, baseline_data){
  
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
        file<- NULL
        print(seed)
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
 
  # totlos.simfs does the same for semi-Markov multistate models using simulation to compute the integral
  # of transition probability since it is not analytically possible
  # returns the expected total time spent in each state starting from the specified state (so it's a vector)
  # I am not able to fix in the simulation scheme the variable that indicates at what age the patient
  # enters the study so, I integrated form 0->t_start and subtracted it to the integration 0->105
  # in Markovian case I forced integration t_start->105 since transition probabilities for t<t_start 
  # where unknown , but I think this is not a bid deal for Semi-Markov case in which the probabilites and
  # integrals are computed from simulation, so data are available for each t>0 because simply simulated
  # from the fitted model
  # (sim.fmsm relies on the presence of a function to sample random numbers from the parametric 
  # survival distribution used in the fitted model x)


  
  # ============================
  # life expectancy ground truth
  # ============================
  
  
  
  cov_means <- colMeans(fits_gompertz_EO[[1]]$data$mml$rate)
  
  newdata <- data.frame(
    cov1 = cov_means[2], 
    cov2 = cov_means[3],
    cov3 = cov_means[4]
  )
  
  tmat <- mstate::transMat(x = list(c(2, 3),c(3),c()), names = c("Dementia-free","Dementia", "Death")) 
 
  gt_tls <-  (totlos.simfs(fits_wei, trans=tmat, newdata = newdata, t=105, cores = cores_lfe) - totlos.simfs(fits_wei, trans=tmat, newdata = newdata, t=t_start, cores = cores_lfe))[1:2]

  #totlos.simfs.mine(fits_wei, tstart= t_start, trans=tmat, newdata = newdata, t=105, cores = cores_lfe)
 
   # ===============
  # EO dataset
  # ===============
  
  flexsurv_tls_EO <- (totlos.simfs(fits_gompertz_EO, trans=tmat, newdata = newdata, t=105, cores = cores_lfe) - totlos.simfs(fits_gompertz_EO, trans=tmat, newdata = newdata, t=t_start, cores = cores_lfe))[1:2]
  bias_flexsurv_tls_EO <- (flexsurv_tls_EO-gt_tls)/gt_tls
  
  # ============
  # coxph
  # ============
  
  if (convergence$coxph[seed]==2) {
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
  
    #basehaz(model_cox,center=FALSE) I get different estimates
    msfit_obj <- msfit(model_cox, newdata = newdata_cox, variance=T, trans=tmat) 
    # confused about times 
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
  flexsurv_tls <-   (totlos.simfs(fits_gompertz, trans=tmat, newdata = newdata, t=105, cores = cores_lfe) - totlos.simfs(fits_gompertz, trans=tmat, newdata = newdata, t=t_start, cores = cores_lfe))[1:2]
  bias_flexsurv_tls <- (flexsurv_tls-gt_tls)/gt_tls
  } else {
    flexsurv_tls <- rep(NA,2)
    bias_flexsurv_tls <- rep(NA,2)
  }
  

  # ======
  # nhm
  # ======
  
  # here we ignore that data comes from a semi-markovian generation process and measure how biased 
  # lfe estimations are once we fit nhm model assuming markovianity 
  
  
  if(convergence$nhm[seed]==2){
    tcrit <- model_nhm$tcrit
    time <- seq(t_start,tcrit-1,by=0.1)
    nhm_probabilities <- predict(model_nhm, times= time)$probabilities 
    nhm_tls<- numeric(ncol(nhm_probabilities))
    
    diff_time <- diff(time)
    nhm_probabilities <- nhm_probabilities[1:length(diff_time),]
    
    for (i in 1:ncol(nhm_probabilities)) {
      nhm_tls[i] <- sum(nhm_probabilities[, i] * diff_time)
    }
    nhm_tls <- nhm_tls[1:2]
    bias_nhm_tls <- (nhm_tls-gt_tls)/gt_tls
  }else{
    nhm_tls <- rep(NA,2)
    bias_nhm_tls <- rep(NA,2)
  }
  


  # ==========
  # imputation
  # ==========
  
  if(convergence$imputation[seed]==2){
    imputation_tls <- c(0,0)
    models_imp <- results_imp[[2]]
    m <- length(models_imp)
    
    for (i in 1:m){
      temp <-  (totlos.simfs(models_imp[[i]], trans=tmat, newdata = newdata, t=105, cores = cores_lfe) - totlos.simfs(models_imp[[i]], trans=tmat, newdata = newdata, t=t_start, cores = cores_lfe))[1:2]
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
    cbind(lfe = nhm_tls, model = "nhm", seed = seed),
    cbind(lfe = imputation_tls, model = "imputation", seed = seed) 
  )
  
  lfe_bias <- rbind(
    cbind(lfe = bias_flexsurv_tls_EO, model = "flexsurv_EO", seed=seed),
    cbind(lfe = bias_coxph_tls, model = "coxph", seed = seed),
    cbind(lfe = bias_flexsurv_tls, model = "flexsurv", seed = seed),
    cbind(lfe = bias_nhm_tls, model = "nhm", seed = seed),
    cbind(lfe = bias_imputation_tls, model = "imputation", seed = seed) 
  )
  
  print(seed)
  return(list(lfe_estimates, lfe_bias))
  
}
 
