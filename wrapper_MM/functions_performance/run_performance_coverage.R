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
  
  min_age <- min(model.msm_age$data[[1]]$age)   
  max_age <- max(model.msm_age$data[[1]]$age)   
  haz <- hazards_mine(model.msm_age, b.covariates = list(age = 0, cov1 = 0, cov2 = 0, cov3 = 0), no.years = 40, CI = F)
  
  for (i in 1:3){
    haz_estimates[i] <- haz[[i]]
  }
  
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
  
  col1 <- numeric()
  col2 <- numeric()
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
  
  # coverage_nhm <- matrix(0, nrow = 3, ncol = 5)
  # 
  # coef_estimates <- model_nhm$par 
  # std_errors <- sqrt(diag(vcov(model_nhm))) 
  # 
  # # Calculate the confidence intervals (95% CI)
  # ci_lower <- coef_estimates - 1.96 * std_errors
  # ci_upper <- coef_estimates + 1.96 * std_errors
  # 
  # # Combine into a data frame for better readability
  # ci <- data.frame(
  #   Estimate = coef_estimates,
  #   Lower = ci_lower,
  #   Upper = ci_upper
  # )
  # params_nhm <- matrix(model_nhm$par, nrow = 3, ncol = 5) 
  # params_nhm <- cbind(params_nhm, exp(params_nhm[,3]), exp(params_nhm[,4]), exp(params_nhm[,5]))
  # colnames(params_nhm) <- colnames(ground_truth_params)
  # 
  # bias_nhm <- compute_bias(params_nhm, ground_truth_params)
  
  # ============
  # imputation
  # ============
  
  coverage_imputation <- matrix(0, nrow = 3, ncol = 5)
  avg_parameters <- results_imp[[1]]
  all_fits <- results_imp[[2]]
  m <- length(all_fits)
  
  ## Let's compute parameters covariance with Rubin rule
  param_matrix <- list()
  for (i in 1:m){
    param_matrix[[i]] <- sapply(all_fits[[i]], coefficients) 
  }
  
  # mean within variance (mean of covariance matrix of each imputation)
  U_list <- list()
  for (i in 1:m){
    U_list[[i]] <- lapply(all_fits[[i]], function(fit) vcov(fit))  
  }
  
  U_bar <- list()
  for (i in 1:3){
    U_bar[[i]] <- matrix(0, nrow = nrow(U_list[[1]][[i]]), ncol = ncol(U_list[[1]][[i]]))
    
    for (j in 1:m){
      U_bar[[i]] <- U_bar[[i]]+U_list[[j]][[i]]
    }
    U_bar[[i]] <- U_bar[[i]]/m
  }
  
  # between variance (variance of parameters estimate)
  param_matrix_mod1 <- list()
  param_matrix_mod2 <- list()
  param_matrix_mod3 <- list()
  for (i in 1:m){
    param_matrix_mod1[[i]] <- param_matrix[[i]][,1]
    param_matrix_mod2[[i]] <- param_matrix[[i]][,2]
    param_matrix_mod3[[i]] <- param_matrix[[i]][,3]
  }
  
  B1 <- rep(0, 5, 5)
  for (j in 1:m) {
    B1 <- B1 + (param_matrix_mod1[[j]] - avg_parameters[1,])%*% t(param_matrix_mod1[[j]] - avg_parameters[1,])
  }
  B1 <- B1 / (m - 1)
  
  B2 <- rep(0, 5, 5)
  for (j in 1:m) {
    B2 <- B2 + (param_matrix_mod1[[j]] - avg_parameters[2,])%*% t(param_matrix_mod1[[j]] - avg_parameters[2,])
  }
  B2 <- B2 / (m - 1)
  
  B3 <- rep(0, 5 ,5)
  for (j in 1:m) {
    B3 <- B3 + (param_matrix_mod1[[j]] - avg_parameters[3,])%*% t(param_matrix_mod1[[j]] - avg_parameters[3,])
  }
  B3 <- B3 / (m - 1)
  
  
  #Rubin's rule
  T1 <- U_bar[[1]] + (1 + 1/m) * B1
  T2 <- U_bar[[2]] + (1 + 1/m) * B2
  T3 <- U_bar[[3]] + (1 + 1/m) * B3
  
  # confidence interval computation
  
  CI_mod1<- compute_CI(T1,U_bar[[1]],B1,m,avg_parameters[1,])
  CI_mod2<- compute_CI(T2,U_bar[[2]],B2,m,avg_parameters[2,])
  CI_mod3<- compute_CI(T3,U_bar[[3]],B3,m,avg_parameters[3,])
  CI_mod1 <- CI_mod1[c(2,1,3,4,5),]
  CI_mod2 <- CI_mod2[c(2,1,3,4,5),]
  CI_mod3 <- CI_mod3[c(2,1,3,4,5),]
  
  for (j in 1:5) {
    ci_lower <- CI_mod1[j, 1]
    ci_upper <- CI_mod1[j, 2]
    coverage_imputation[1,j] <- (ground_truth_params[1,j] >= ci_lower) && (ground_truth_params[1,j] <= ci_upper)
  }
  for (j in 1:5) {
    ci_lower <- CI_mod2[j, 1]
    ci_upper <- CI_mod2[j, 2]
    coverage_imputation[2,j] <- (ground_truth_params[2,j] >= ci_lower) && (ground_truth_params[2,j] <= ci_upper)
  }
  for (j in 1:5) {
    ci_lower <- CI_mod3[j, 1]
    ci_upper <- CI_mod3[j, 2]
    coverage_imputation[3,j] <- (ground_truth_params[3,j] >= ci_lower) && (ground_truth_params[3,j] <= ci_upper)
  }
  
  param_names <- colnames(ground_truth_params)
  colnames(coverage_imputation) <- param_names
  
  
  coverage_tot <- rbind(
    cbind(coverage_coxph, model = "coxph", seed = seed, transition = c(1,2,3)),
    cbind(coverage_flexsurv, model = "flexusrv", seed = seed, transition = c(1,2,3)),
    cbind(coverage_msm, model = "msm", seed = seed, transition = c(1,2,3)),
    cbind(coverage_msm_age, model = "msm_age", seed = seed, transition = c(1,2,3)),
    #cbind(coverage_nhm, model = "nhm", seed = seed, transition = c(1,2,3)),
    cbind(coverage_imputation, model = "imputation", seed = seed, transition = c(1,2,3))
  )

  return (bias_tot)
}

