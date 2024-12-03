computing_life_expectancy <- function(n_pats, scheme, seed, convergence, t_start, covs){
  
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
                       #"model_nhm.RData", 
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
  
  #####################
  # trying with hesim
  #####################
  
  #passare meanlog e sdlog da esterno come per covs
  meanlog <- mean(log(fits_wei[[1]]$data$mml$rate[,2]))
  sdlog <- sd(log(fits_wei[[1]]$data$mml$rate[,2]))
  
  #for each model i should pass different fits Idk about cox
  #what about censoring?
  sim_data <- simulation(n_pats, fits, meanlog, sdlog, covs)
  
  temp <- dismod$sim_disease(max_t=120, max_age=120)
  # I get prob of being in each state in each age how to compute lfe??
  t <- seq(60,105,0.1)
  probs <- dismod$sim_stateprobs(t)
  probs[, c("sample", "strategy_id", "grp_id") := NULL]
  
  # ============================
  # life expectancy ground truth
  # ============================
  
  
  tmat <- mstate::transMat(x = list(c(2, 3),c(3),c()), names = c("Dementia-free","Dementia", "Death")) 
 
  gt_tls <-  (totlos.simfs(fits_wei, trans=tmat, newdata = covs, t=105, cores = cores_lfe) - totlos.simfs(fits_wei, trans=tmat, newdata = newdata, t=t_start, cores = cores_lfe))[1:2]

  #totlos.simfs.mine(fits_wei, tstart= t_start, trans=tmat, newdata = newdata, t=105, cores = cores_lfe)
 
   # ===============
  # EO dataset
  # ===============
  
  flexsurv_tls_EO <- (totlos.simfs(fits_gompertz_EO, trans=tmat, newdata = covs, t=105, cores = cores_lfe) - totlos.simfs(fits_gompertz_EO, trans=tmat, newdata = newdata, t=t_start, cores = cores_lfe))[1:2]
  bias_flexsurv_tls_EO <- (flexsurv_tls_EO-gt_tls)/gt_tls
  
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
  
    #basehaz(model_cox,center=FALSE) I get different estimates
    msfit_obj <- msfit(model_cox, newdata = newdata_cox, variance=T, trans=tmat) 
    # confused about times
    tv <- unique(fit$Haz$time)
    #forse non devo usare reset per come ho preparato i dati
    cox_trans_prob <- mssample(msfit_obj$Haz, tmat, history = list(state = 1, time = min(msfit_obj$time)), )
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
  flexsurv_tls <-   (totlos.simfs(fits_gompertz, trans=tmat, newdata = covs, t=105, cores = cores_lfe) - totlos.simfs(fits_gompertz, trans=tmat, newdata = newdata, t=t_start, cores = cores_lfe))[1:2]
  bias_flexsurv_tls <- (flexsurv_tls-gt_tls)/gt_tls
  } else {
    flexsurv_tls <- rep(NA,2)
    bias_flexsurv_tls <- rep(NA,2)
  }
  


  # ==========
  # imputation
  # ==========
  
  if(convergence$imputation[seed]==2){
    imputation_tls <- c(0,0)
    models_imp <- results_imp[[2]]
    m <- length(models_imp)
    
    for (i in 1:m){
      temp <-  (totlos.simfs(models_imp[[i]], trans=tmat, newdata = covs, t=105, cores = cores_lfe) - totlos.simfs(models_imp[[i]], trans=tmat, newdata = newdata, t=t_start, cores = cores_lfe))[1:2]
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
 
