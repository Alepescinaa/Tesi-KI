width_ic <- function(n_pats, scheme, seed, convergence) {

  
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
                       "results_imp.RData",
                       "computational_time.RData")
    
    for (file in files_to_load) {
      if (file.exists(file)) {
        load(file)
      } else {
        warning(paste("File does not exist:", file))
        file <- sub("\\.Rdata$", "", file, ignore.case = T) 
        file<- NULL
      }
    }
  } else {
    warning(paste("Seed directory does not exist:", seed_dir))
  }
  
  # ============
  # EO dataset
  # ============
  
  width_EO <- matrix(0, nrow = 3, ncol = 3)
  
  ci_lower <- lapply(fits_gompertz_EO, function(model) confint(model)[, 1])
  ci_upper <- lapply(fits_gompertz_EO, function(model) confint(model)[, 2])
  
  ci_lower_EO <- do.call(cbind, ci_lower)
  ci_upper_EO <- do.call(cbind, ci_upper)
  ci_lower_EO <- ci_lower_EO[c(3, 4, 5), ]
  ci_upper_EO <- ci_upper_EO[c(3, 4, 5), ]
  
  width_EO <- round(ci_upper_EO-ci_lower_EO,3)
  rownames(width_EO) <- c("cov1","cov2","cov3")
  colnames(width_EO) <- c("1","2", "3")
  
  rownames(ci_upper_EO) <- c("cov1","cov2","cov3")
  colnames(ci_upper_EO) <- c("1","2", "3")
  rownames(ci_lower_EO) <- c("cov1","cov2","cov3")
  colnames(ci_lower_EO) <- c("1","2", "3")
  
  # =========
  # coxph
  # =========
  
  if (convergence$coxph[seed]==2){
    width_cox <- matrix(0, nrow = 3, ncol = 3)
    
    ci_lower <- confint(model_cox)[,1]
    ci_upper <- confint(model_cox)[,2]
    
    ci_lower_cox <- matrix(ci_lower,3,3, byrow=F)
    ci_upper_cox <- matrix(ci_upper,3,3, byrow=F)
    
    width_cox <- round(ci_upper_cox-ci_lower_cox,3)
    rownames(width_cox) <- c("cov1", "cov2", "cov3")
    colnames(width_cox) <- c("1", "2", "3")
    
    rownames(ci_upper_cox) <- c("cov1", "cov2", "cov3")
    colnames(ci_upper_cox) <- c("1", "2", "3")
    rownames(ci_lower_cox) <- c("cov1", "cov2", "cov3")
    colnames(ci_lower_cox) <- c("1", "2", "3")
    
  } else{
    width_cox <- matrix(NA, nrow = 3, ncol = 3)
    rownames(width_cox) <- c("cov1", "cov2", "cov3")
    colnames(width_cox) <- c("1", "2", "3")  
    
    ci_upper_cox <- matrix(NA, nrow = 3, ncol = 3)
    rownames(ci_upper_cox) <- c("cov1", "cov2", "cov3")
    colnames(ci_upper_cox) <- c("1", "2", "3")  
    
    ci_lower_cox <- matrix(NA, nrow = 3, ncol = 3)
    rownames(ci_lower_cox) <- c("cov1", "cov2", "cov3")
    colnames(ci_lower_cox) <- c("1", "2", "3")  
    }
  
  
  # ============ 
  # flexsurv 
  # ============
  
  if(convergence$flexsurv[seed]==2){
    width_flexsurv <- matrix(0, nrow = 3, ncol = 3)
    
    ci_lower <- lapply(fits_gompertz, function(model) confint(model)[, 1])
    ci_upper <- lapply(fits_gompertz, function(model) confint(model)[, 2])
    
    ci_lower_flex <- do.call(cbind, ci_lower)
    ci_upper_flex <- do.call(cbind, ci_upper)
    ci_lower_flex <- ci_lower_flex[c(3, 4, 5), ]
    ci_upper_flex <- ci_upper_flex[c(3, 4, 5), ]
    
    width_flexsurv <- round(ci_upper_flex-ci_lower_flex,3)
    rownames(width_flexsurv) <- c("cov1", "cov2", "cov3")
    colnames(width_flexsurv) <- c("1", "2", "3")  
    
    rownames(ci_upper_flex) <- c("cov1", "cov2", "cov3")
    colnames(ci_upper_flex) <- c("1", "2", "3")
    rownames(ci_lower_flex) <- c("cov1", "cov2", "cov3")
    colnames(ci_lower_flex) <- c("1", "2", "3")
    
  } else{
    width_flexsurv <- matrix(NA, nrow = 3, ncol = 3)
    rownames(width_flexsurv) <- c("cov1", "cov2", "cov3")
    colnames(width_flexsurv) <- c("1", "2", "3")
    
    ci_upper_flex <- matrix(NA, nrow = 3, ncol = 3)
    rownames(ci_upper_flex) <- c("cov1", "cov2", "cov3")
    colnames(ci_upper_flex) <- c("1", "2", "3")
    
    ci_lower_flex <- matrix(NA, nrow = 3, ncol = 3)
    rownames(ci_lower_flex) <- c("cov1", "cov2", "cov3")
    colnames(ci_lower_flex) <- c("1", "2", "3")
  }
  
  
  # ============ 
  # imputation 
  # ============ 
  
  if(convergence$imputation[seed]==2){
    width_imputation <- matrix(0, nrow = 3, ncol = 3)
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
    
    CI_mod1 <- CI_mod1[c( 3, 4, 5), ]
    CI_mod2 <- CI_mod2[c( 3, 4, 5), ]
    CI_mod3 <- CI_mod3[c( 3, 4, 5), ]
    
    ci_lower_imp <- cbind(CI_mod1[, 1], CI_mod2[, 1], CI_mod3[, 1])
    ci_upper_imp <- cbind(CI_mod1[, 2], CI_mod2[, 2], CI_mod3[, 2])
    
    width_imputation <- round(ci_upper_imp-ci_lower_imp,3)
    rownames(width_imputation) <- c("cov1", "cov2", "cov3")
    colnames(width_imputation) <- c("1", "2", "3")  
    
    rownames(ci_upper_imp) <- c("cov1", "cov2", "cov3")
    colnames(ci_upper_imp) <- c("1", "2", "3")  
    rownames(ci_lower_imp) <- c("cov1", "cov2", "cov3")
    colnames(ci_lower_imp) <- c("1", "2", "3")  
    
  } else {
    width_imputation <- matrix(NA, nrow = 3, ncol = 3)
    rownames(width_imputation) <- c("cov1", "cov2", "cov3")
    colnames(width_imputation) <- c("1", "2", "3")  
    
    ci_upper_imp <- matrix(NA, nrow = 3, ncol = 3)
    rownames(ci_upper_imp) <- c("cov1", "cov2", "cov3")
    colnames(ci_upper_imp) <- c("1", "2", "3")  
    
    ci_lower_imp <- matrix(NA, nrow = 3, ncol = 3)
    rownames(ci_lower_imp) <- c("cov1", "cov2", "cov3")
    colnames(ci_lower_imp) <- c("1", "2", "3")  
  }
  
  width_tot <- rbind(
    cbind(width_EO, model = "flexsurv_EO", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(width_cox, model = "coxph", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(width_flexsurv, model = "flexsurv", seed = seed, covariate=c("cov1","cov2","cov3")),
    #cbind(width_nhm, model = "nhm", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(width_imputation, model = "imputation", seed = seed, covariate=c("cov1","cov2","cov3"))
  )
  
  lower_tot <- rbind(
    cbind(ci_lower_EO, model = "flexsurv_EO", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(ci_lower_cox, model = "coxph", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(ci_lower_flex, model = "flexsurv", seed = seed, covariate=c("cov1","cov2","cov3")),
    #cbind(ci_lower_nhm, model = "nhm", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(ci_lower_imp, model = "imputation", seed = seed, covariate=c("cov1","cov2","cov3"))
  )
  
  upper_tot <- rbind(
    cbind(ci_upper_EO, model = "flexsurv_EO", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(ci_upper_cox, model = "coxph", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(ci_upper_flex, model = "flexsurv", seed = seed, covariate=c("cov1","cov2","cov3")),
   # cbind(ci_upper_nhm, model = "nhm", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(ci_upper_imp, model = "imputation", seed = seed, covariate=c("cov1","cov2","cov3"))
  )
  
  lower_tot <- as.data.frame(lower_tot)
  upper_tot <- as.data.frame(upper_tot)
  lower_tot <- lower_tot %>% mutate(ic = "lower")
  upper_tot <- upper_tot %>% mutate(ic = "upper")
  ic <- full_join(lower_tot,upper_tot)
  
  return(ic)
}

