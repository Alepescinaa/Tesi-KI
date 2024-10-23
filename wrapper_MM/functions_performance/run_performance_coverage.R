run_performance_coverage <- function(n_pats, scheme, seed){
  
  load("ground_truthMM.RData")
  ground_truth_params <- ground_truth_params[,1:5]
  
  if (n_pats == 500){
    if (scheme == 2){
      scheme_dir <- "results_500/saved_models_scheme2"
    } else if (scheme == 3){
      scheme_dir <-  "results_500/saved_models_scheme3"
    } else if (scheme == 4){
      scheme_dir <-  "results_500/saved_models_scheme4"
    } else if (scheme == 5){
      scheme_dir <-  "results_500/saved_models_scheme5"
    }
  } else if (n_pats == 2000){
    if (scheme == 2){
      scheme_dir <- "results_2K/saved_models_scheme2"
    } else if (scheme == 3){
      scheme_dir <-  "results_2K/saved_models_scheme3"
    } else if (scheme == 4){
      scheme_dir <-  "results_2K/saved_models_scheme4"
    } else if (scheme == 5){
      scheme_dir <-  "results_2K/saved_models_scheme5"
    }
  } else if (n_pats == 5000){
    if (scheme == 2){
      scheme_dir <- "results_5K/saved_models_scheme2"
    } else if (scheme == 3){
      scheme_dir <-  "results_5K/saved_models_scheme3"
    } else if (scheme == 4){
      scheme_dir <-  "results_5K/saved_models_scheme4"
    } else if (scheme == 5){
      scheme_dir <-  "results_5K/saved_models_scheme5"
    }
  } else if (n_pats == 10000){
    if (scheme == 2){
      scheme_dir <- "results_10K/saved_models_scheme2"
    } else if (scheme == 3){
      scheme_dir <-  "results_10K/saved_models_scheme3"
    } else if (scheme == 4){
      scheme_dir <-  "results_10K/saved_models_scheme4"
    } else if (scheme == 5){
      scheme_dir <-  "results_10K/saved_models_scheme5"
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
      }
    }
  } else {
    warning(paste("Seed directory does not exist:", seed_dir))
  }
  
  
  # =========
  # coxph
  # =========
  
  coverage_cox <- matrix(0, nrow = 3, ncol = 3)
  
  
  for (i in 1:3){
      ci<- confint(model_cox[[i]])
      for (j in 1:length(ci[,1])) {
        ci_lower <- ci[j, 1]
        ci_upper <- ci[j, 2]
        coverage_cox[i,j] <- (ground_truth_params[i,j+2] >= ci_lower) && (ground_truth_params[i,j+2] <= ci_upper)
      }
  }
  
  empty <- numeric(3)
  empty[] <- NA  
  coverage_cox <- cbind(empty, empty,coverage_cox)
  param_names <- colnames(ground_truth_params)
  colnames(coverage_cox) <- param_names
  
  # ============
  # flexsurv
  # ============
  
  coverage_flexsurv <- matrix(0, nrow = 3, ncol = 5)
 
  for (i in 1:3){
    for (j in 1:5){
      ci<- confint(fits_gompertz[[i]])
      ci <- ci[c(2,1,3,4,5),]
      for (j in 1:length(ci[,1])) {
        ci_lower <- ci[j, 1]
        ci_upper <- ci[j, 2]
        coverage_flexsurv[i,j] <- (ground_truth_params[i,j] >= ci_lower) && (ground_truth_params[i,j] <= ci_upper)
      }
    }
  } 
 
  param_names <- colnames(ground_truth_params)
  colnames(coverage_flexsurv) <- param_names
  
  # ============
  # msm
  # ============
  
  coverage_msm <- matrix(0, nrow = 3, ncol = 3)
  
  for (k in 1:3) {
    ci <- model.msm$ci[(3 * (k - 1) + 4):(3 * k + 3), ] 
    for (j in 1:nrow(ci)) {
      ci_lower <- ci[j, 1]
      ci_upper <- ci[j, 2]
      coverage_msm[k, j] <- (ground_truth_params[k, j + 2] >= ci_lower) && (ground_truth_params[k, j + 2] <= ci_upper)
    }
  }
  
  empty <- numeric(3)
  empty[] <- NA  
  coverage_msm <- cbind(empty, empty,coverage_msm)
  param_names <- colnames(ground_truth_params)
  colnames(coverage_msm) <- param_names 
  

  # ============
  # msm + age
  # ============
  
  # FIX RATE E SHAPE CI
  

  min_age <- min(model.msm_age$data[[1]]$age)   
  max_age <- max(model.msm_age$data[[1]]$age)   
  haz <- hazards_mine(model.msm_age, b.covariates = list(age = 0, cov1 = 0, cov2 = 0, cov3 = 0), no.years = 40, CI = T)
  
  haz_estimates <- list()
  for (i in 1:3){
    haz_estimates[[i]] <- haz[[1]][[i]]
  }
 
  # rate <- numeric(3)   
  # shape <- numeric(3)   
  # rate_CI <- matrix(NA, nrow = 3, ncol = 2) 
  # shape_CI <- matrix(NA, nrow = 3, ncol = 2)
  # 
  # for (i in 1:3) {   
  #   age_grid <- seq(min_age, max_age, length.out = length(unlist(haz_estimates[[1]])))   
  #   y <- log(as.numeric(unlist(haz_estimates[[i]])))  
  #   reg_model <- lm(y ~ age_grid)   
  #   rate[i] <- reg_model$coefficients[1]   
  #   shape[i] <- reg_model$coefficients[2] 
  # }
  # for (i in 1:3) {   
  #   y <- log(as.numeric(unlist(haz_LB[[i]])))  
  #   reg_model <- lm(y ~ age_grid)   
  #   rate_CI[i,1] <- reg_model$coefficients[1]   
  #   shape_CI[i,1] <- reg_model$coefficients[2] 
  # }
  # for (i in 1:3) {   
  #   y <- log(as.numeric(unlist(haz_UB[[i]])))  
  #   reg_model <- lm(y ~ age_grid)   
  #   rate_CI[i,2] <- reg_model$coefficients[1]   
  #   shape_CI[i,2] <- reg_model$coefficients[2] 
  # }
  #   
  
  set.seed(2024) 
  n_bootstrap <- 1000 
  rate_bootstrap <- matrix(NA, nrow = n_bootstrap, ncol = 3)
  shape_bootstrap <- matrix(NA, nrow = n_bootstrap, ncol = 3)
  
  for (i in 1:3) {
    for (b in 1:n_bootstrap) {
    
      sample_indices <- sample(1:length(haz_estimates[[i]]), replace = TRUE)
      sample_haz <- unlist(haz_estimates[[i]])[sample_indices]
      
      age_grid <- seq(min_age, max_age, length.out = length(sample_haz))
      y <- log(sample_haz)
      reg_model <- lm(y ~ age_grid)
      
      rate_bootstrap[b, i] <- reg_model$coefficients[1]
      shape_bootstrap[b, i] <- reg_model$coefficients[2]
    }
  }
  
  rate_CI <- apply(rate_bootstrap, 2, quantile, probs = c(0.025, 0.975))
  shape_CI <- apply(shape_bootstrap, 2, quantile, probs = c(0.025, 0.975))
  
  coverage_msm_age <- matrix(0, nrow = 3, ncol = 3)
  
  ci <- model.msm_age$ci[4:6, ] 
  for (j in 1:nrow(ci)) {
      ci_lower <- ci[j, 1]
      ci_upper <- ci[j, 2]
      coverage_msm_age[1, j] <- (ground_truth_params[1, j + 2] >= ci_lower) && (ground_truth_params[1, j + 2] <= ci_upper)
    }
  ci <- model.msm_age$ci[8:10, ] 
  for (j in 1:nrow(ci)) {
  ci_lower <- ci[j, 1]
  ci_upper <- ci[j, 2]
  coverage_msm_age[2, j] <- (ground_truth_params[2, j + 2] >= ci_lower) && (ground_truth_params[2, j + 2] <= ci_upper)
  }
  ci <- model.msm_age$ci[12:14, ] 
  for (j in 1:nrow(ci)) {
  ci_lower <- ci[j, 1]
  ci_upper <- ci[j, 2]
  coverage_msm_age[3, j] <- (ground_truth_params[3, j + 2] >= ci_lower) && (ground_truth_params[3, j + 2] <= ci_upper)
  }
  
  for (i in 1:3){
    col1[i] <- (ground_truth_params[i,1] >= rate_CI[1,i]) && (ground_truth_params[i,1] <= rate_CI[2,i])
    col2[i] <- (ground_truth_params[i,2] >= shape_CI[1,i]) && (ground_truth_params[i,2] <= shape_CI[2,i])
  }

  coverage_msm_age <- cbind(col1,col2,coverage_msm_age)
  param_names <- colnames(ground_truth_params)
  colnames(coverage_msm_age) <- param_names 
  
  # ============
  # nhm
  # ============
  
  coverage_nhm <- matrix(0, nrow = 3, ncol = 5)
  
  coef_estimates <- model_nhm$par 
  std_errors <- sqrt(diag(vcov(model_nhm))) 
  
  # Calculate the confidence intervals (95% CI)
  ci_lower <- coef_estimates - 1.96 * std_errors
  ci_upper <- coef_estimates + 1.96 * std_errors
  
  # Combine into a data frame for better readability
  ci <- data.frame(
    Estimate = coef_estimates,
    Lower = ci_lower,
    Upper = ci_upper
  )
  params_nhm <- matrix(model_nhm$par, nrow = 3, ncol = 5) 
  params_nhm <- cbind(params_nhm, exp(params_nhm[,3]), exp(params_nhm[,4]), exp(params_nhm[,5]))
  colnames(params_nhm) <- colnames(ground_truth_params)
  
  bias_nhm <- compute_bias(params_nhm, ground_truth_params)
  
  # ============
  # imputation
  # ============
  
  coverage_imputation <- matrix(0, nrow = 3, ncol = 5)
  results_imp
  for (i in 1:3){
    for (j in 1:5){
      ci<- confint(fits_gompertz[[i]])
      ci <- ci[c(2,1,3,4,5),]
      for (j in 1:length(ci[,1])) {
        ci_lower <- ci[j, 1]
        ci_upper <- ci[j, 2]
        coverage_imputation[i,j] <- (ground_truth_params[i,j] >= ci_lower) && (ground_truth_params[i,j] <= ci_upper)
      }
    }
  } 
  
  param_names <- colnames(ground_truth_params)
  colnames(coverage_imputation) <- param_names
  
  
  # bias_tot <- rbind(
  #   cbind(bias_coxph, model = "coxph", seed = seed, transition = c(1,2,3)),
  #   cbind(bias_flexsurv, model = "flexusrv", seed = seed, transition = c(1,2,3)),
  #   cbind(bias_msm, model = "msm", seed = seed, transition = c(1,2,3)),
  #   cbind(bias_msm_age, model = "msm_age", seed = seed, transition = c(1,2,3)),
  #   cbind(bias_nhm, model = "nhm", seed = seed, transition = c(1,2,3)),
  #   cbind(bias_imputation, model = "imputation", seed = seed, transition = c(1,2,3))
  # )
  # 
  return (bias_tot)
}

