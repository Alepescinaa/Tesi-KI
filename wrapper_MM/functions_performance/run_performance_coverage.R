run_performance_coverage <- function(n_pats, scheme, seed, convergence) {
  load("./wrapper_MM/ground_truthMM.RData")
  ground_truth_params <- ground_truth_params[, 1:5]
  
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
  
  # =========
  # coxph
  # =========
  
  if (convergence$coxph[seed]==2){
    coverage_cox <- matrix(0, nrow = 3, ncol = 5)
    
    ci_lower <- lapply(model_cox, function(model) confint(model)[, 1])
    ci_upper <- lapply(model_cox, function(model) confint(model)[, 2])
    
    ci_lower_cox <- do.call(cbind, ci_lower)
    ci_upper_cox <- do.call(cbind, ci_upper)
    
    coverage_cox <- compute_coverage(ci_lower_cox, ci_upper_cox, ground_truth_params)
    coverage_cox[,1:2] <- NA
  } else{
    coverage_cox <- matrix(NA, nrow = nrow(ground_truth_params), ncol = ncol(ground_truth_params))
    colnames(coverage_cox) <- colnames(ground_truth_params)
    rownames(coverage_cox) <- rownames(ground_truth_params)
  }
  
  
  # ============ 
  # flexsurv 
  # ============
  
  if(convergence$flexsurv[seed]==2){
    coverage_flexsurv <- matrix(0, nrow = 3, ncol = 5)
    
    ci_lower <- lapply(fits_gompertz, function(model) confint(model)[, 1])
    ci_upper <- lapply(fits_gompertz, function(model) confint(model)[, 2])
    
    ci_lower_flex <- do.call(cbind, ci_lower)
    ci_upper_flex <- do.call(cbind, ci_upper)
    ci_lower_flex <- ci_lower_flex[c(2, 1, 3, 4, 5), ]
    ci_upper_flex <- ci_upper_flex[c(2, 1, 3, 4, 5), ]
    
    coverage_flexsurv <- compute_coverage(ci_lower_flex, ci_upper_flex, ground_truth_params)
  } else{
    coverage_flexsurv <- matrix(NA, nrow = nrow(ground_truth_params), ncol = ncol(ground_truth_params))
    colnames(coverage_flexsurv) <- colnames(ground_truth_params)
    rownames(coverage_flexsurv) <- rownames(ground_truth_params)
  }
  
  
  # ============
  # msm
  # ============
  
  if (convergence$msm[seed]==2){
    coverage_msm <- matrix(0, nrow = 3, ncol = 5)
    
    ci_lower_msm <- model.msm$ci[4:12, 1] # are ordered by transition ex 1.cov1,2.cov1,3.cov1,1.cov2,2.cov2..
    ci_upper_msm <- model.msm$ci[4:12, 2]
    
    ci_lower_msm <- matrix(ci_lower_msm, nrow = 3, ncol = 3, byrow= T)
    ci_upper_msm <- matrix(ci_upper_msm, nrow = 3, ncol = 3, byrow= T)
    rownames(ci_lower_msm) <- c("cov1", "cov2", "cov3")
    rownames(ci_upper_msm) <- c("cov1", "cov2", "cov3")
    
    coverage_msm <- compute_coverage(ci_lower_msm, ci_upper_msm, ground_truth_params)
    coverage_msm[,1:2] <- NA
  } else{
    coverage_msm <- matrix(NA, nrow = nrow(ground_truth_params), ncol = ncol(ground_truth_params))
    colnames(coverage_msm) <- colnames(ground_truth_params)
    rownames(coverage_msm) <- rownames(ground_truth_params)
  }
 
  
  # ============
  # msm + age
  # ============
  
  if (convergence$msm_age[seed]==2){
    min_age <- min(model.msm_age$data[[1]]$age)
    max_age <- max(model.msm_age$data[[1]]$age)
    haz <- hazards_mine(model.msm_age, b.covariates = list(age = 0, cov1 = 0, cov2 = 0, cov3 = 0), no.years = 40, CI = F)
    
    haz_estimates <- numeric()
    for (i in 1:3) {
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
    
    ci_lower_msm_age <- model.msm_age$ci[c(4:6, 8:10, 12:14), 1]
    ci_upper_msm_age <- model.msm_age$ci[c(4:6, 8:10, 12:14), 2]
    
    ci_lower_msm_age <- matrix(ci_lower_msm_age, nrow = 3, ncol = 3, byrow=T)
    ci_upper_msm_age <- matrix(ci_upper_msm_age, nrow = 3, ncol = 3, byrow=T)
    
    ci_lower_msm_age <- rbind(rate_CI[1, ], shape_CI[1, ], ci_lower_msm_age)
    ci_upper_msm_age <- rbind(rate_CI[2, ], shape_CI[2, ], ci_upper_msm_age)
    rownames(ci_lower_msm_age) <- c("rate", "shape", "cov1", "cov2", "cov3")
    rownames(ci_upper_msm_age) <- c("rate", "shape", "cov1", "cov2", "cov3")
    
    coverage_msm_age <- compute_coverage(ci_lower_msm_age, ci_upper_msm_age, ground_truth_params)
  } else{
    coverage_msm_age <- matrix(NA, nrow = nrow(ground_truth_params), ncol = ncol(ground_truth_params))
    colnames(coverage_msm_age) <- colnames(ground_truth_params)
    rownames(coverage_msm_age) <- rownames(ground_truth_params)
  }
  
  
  # ========
  # nhm
  # ========
  
  if (convergence$nhm[seed]==2){
    ci_nhm <- get_params_nhm(model_nhm, ci = TRUE)
    
    ci_lower_nhm <- ci_nhm[,2]
    ci_lower_nhm <- matrix(ci_lower_nhm, 5, 3, byrow = T)
    ci_upper_nhm <- ci_nhm[,3]
    ci_upper_nhm <- matrix(ci_upper_nhm, 5, 3, byrow = T)
    rownames(ci_lower_nhm) <- c("rate", "shape", "cov1", "cov2", "cov3")
    rownames(ci_upper_nhm) <- c("rate", "shape", "cov1", "cov2", "cov3")
    
    coverage_nhm <- compute_coverage(ci_lower_nhm, ci_upper_nhm, ground_truth_params)
  } else {
    coverage_nhm <- matrix(NA, nrow = nrow(ground_truth_params), ncol = ncol(ground_truth_params))
    colnames(coverage_nhm) <- colnames(ground_truth_params)
    rownames(coverage_nhm) <- rownames(ground_truth_params)
    }
  
  # ============ 
  # imputation 
  # ============ 
  
  if(convergence$imputation[seed]==2){
    coverage_imputation <- matrix(0, nrow = 3, ncol = 5)
    avg_parameters <- results_imp[[1]]
    all_fits <- results_imp[[2]]
    m <- length(all_fits)
    
    # Let's compute parameters covariance with Rubin rule
    param_matrix <- list()
    for (i in 1:m) {
      param_matrix[[i]] <- sapply(all_fits[[i]], coefficients)
    }
    
    # mean within variance (mean of covariance matrix of each imputation)
    U_list <- list()
    for (i in 1:m) {
      U_list[[i]] <- lapply(all_fits[[i]], function(fit) vcov(fit))
    }
    
    U_bar <- list()
    for (i in 1:3) {
      U_bar[[i]] <- matrix(0, nrow = nrow(U_list[[1]][[i]]), ncol = ncol(U_list[[1]][[i]]))
      
      for (j in 1:m) {
        U_bar[[i]] <- U_bar[[i]] + U_list[[j]][[i]]
      }
      U_bar[[i]] <- U_bar[[i]] / m
    }
    
    # between variance (variance of parameters estimate)
    param_matrix_mod1 <- list()
    param_matrix_mod2 <- list()
    param_matrix_mod3 <- list()
    for (i in 1:m) {
      param_matrix_mod1[[i]] <- param_matrix[[i]][, 1]
      param_matrix_mod2[[i]] <- param_matrix[[i]][, 2]
      param_matrix_mod3[[i]] <- param_matrix[[i]][, 3]
    }
    
    B1 <- matrix(0, nrow = 5, ncol = 5)
    for (j in 1:m) {
      B1 <- B1 + (param_matrix_mod1[[j]] - avg_parameters[1, ]) %*% t(param_matrix_mod1[[j]] - avg_parameters[1, ])
    }
    B1 <- B1 / (m - 1)
    
    B2 <- matrix(0, nrow = 5, ncol = 5)
    for (j in 1:m) {
      B2 <- B2 + (param_matrix_mod2[[j]] - avg_parameters[2, ]) %*% t(param_matrix_mod2[[j]] - avg_parameters[2, ])
    }
    B2 <- B2 / (m - 1)
    
    B3 <- matrix(0, nrow = 5, ncol = 5)
    for (j in 1:m) {
      B3 <- B3 + (param_matrix_mod3[[j]] - avg_parameters[3, ]) %*% t(param_matrix_mod3[[j]] - avg_parameters[3, ])
    }
    B3 <- B3 / (m - 1)
    
    # Rubin's rule
    T1 <- U_bar[[1]] + (1 + 1/m) * B1
    T2 <- U_bar[[2]] + (1 + 1/m) * B2
    T3 <- U_bar[[3]] + (1 + 1/m) * B3
    
    # confidence interval computation
    CI_mod1 <- compute_CI(T1, U_bar[[1]], B1, m, avg_parameters[1, ])
    CI_mod2 <- compute_CI(T2, U_bar[[2]], B2, m, avg_parameters[2, ])
    CI_mod3 <- compute_CI(T3, U_bar[[3]], B3, m, avg_parameters[3, ])
    
    CI_mod1 <- CI_mod1[c(2, 1, 3, 4, 5), ]
    CI_mod2 <- CI_mod2[c(2, 1, 3, 4, 5), ]
    CI_mod3 <- CI_mod3[c(2, 1, 3, 4, 5), ]
    
    ci_lower_imp <- cbind(CI_mod1[, 1], CI_mod2[, 1], CI_mod3[, 1])
    ci_upper_imp <- cbind(CI_mod1[, 2], CI_mod2[, 2], CI_mod3[, 2])
    
    coverage_imputation <- compute_coverage(ci_lower_imp, ci_upper_imp, ground_truth_params)
  } else {
    coverage_imputation <- matrix(NA, nrow = nrow(ground_truth_params), ncol = ncol(ground_truth_params))
    colnames(coverage_imputation) <- colnames(ground_truth_params)
    rownames(coverage_imputation) <- rownames(ground_truth_params)
  }
  
  coverage_tot <- rbind(
    cbind(coverage_cox, model = "coxph", seed = seed, transition = c(1, 2, 3)),
    cbind(coverage_flexsurv, model = "flexsurv", seed = seed, transition = c(1, 2, 3)),
    cbind(coverage_msm, model = "msm", seed = seed, transition = c(1, 2, 3)),
    cbind(coverage_msm_age, model = "msm_age", seed = seed, transition = c(1, 2, 3)),
    cbind(coverage_nhm, model = "nhm", seed = seed, transition = c(1, 2, 3)),
    cbind(coverage_imputation, model = "imputation", seed = seed, transition = c(1, 2, 3))
  )
  
  return(coverage_tot)
}
