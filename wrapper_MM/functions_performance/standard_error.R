standard_error <- function(n_pats, scheme, seed, convergence) {
  
  
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
        print(file)
        print(seed)
        warning(paste("File does not exist:", file))
        model_nhm <- NULL
      }
    }
  } else {
    warning(paste("Seed directory does not exist:", seed_dir))
  }
  # ============
  # EO dataset
  # ============
  
  sd_EO <- matrix(0, nrow = 3, ncol = 3)
  
  for (i in 1:3){
    sd_EO[,i] <- fits_gompertz_EO[[i]]$res.t[(3:5),4]
  }
  
  rownames(sd_EO) <- c("cov1","cov2","cov3")
  colnames(sd_EO) <- c("1","2", "3")
  
  
  # =========
  # coxph
  # =========
  
  if (convergence$coxph[seed]==2){
    sd_cox <- matrix(sqrt(diag(vcov(model_cox))), nrow = 3, ncol = 3, byrow=F)
  
    rownames(sd_cox) <- c("cov1","cov2","cov3")
    colnames(sd_cox) <- c("1","2", "3")
    
  } else{
    sd_cox <- matrix(NA, nrow = 3, ncol = 3)
    rownames(sd_cox) <- c("cov1", "cov2", "cov3")
    colnames(sd_cox) <- c("1", "2", "3")  

  }
  
  
  # ============ 
  # flexsurv 
  # ============
  
  if(convergence$flexsurv[seed]==2){
    sd_flexsurv <- matrix(0, nrow = 3, ncol = 3)
    
    for (i in 1:3){
      sd_flexsurv[,i] <- fits_gompertz[[i]]$res.t[(3:5),4]
    }
    
    rownames(sd_flexsurv) <- c("cov1","cov2","cov3")
    colnames(sd_flexsurv) <- c("1","2", "3")
    
  } else{
    sd_flexsurv <- matrix(NA, nrow = 3, ncol = 3)
    rownames(sd_flexsurv) <- c("cov1", "cov2", "cov3")
    colnames(sd_flexsurv) <- c("1", "2", "3")
  
  }
  
  
  # ============
  # msm
  # ============
  
  if (convergence$msm[seed]==2){
    sd <- sqrt(diag(model.msm$covmat))[4:12]
    sd_msm <- matrix(sd, nrow = 3, ncol = 3, byrow=T)# are ordered by transition ex 1.cov1,2.cov1,3.cov1,1.cov2,2.cov2..
    
    rownames(sd_msm) <- c("cov1", "cov2", "cov3")
    colnames(sd_msm) <- c("1", "2", "3")  
    
  } else{
    sd_msm <- matrix(NA, nrow = 3, ncol = 3)
    rownames(sd_msm) <- c("cov1", "cov2", "cov3")
    colnames(sd_msm) <- c("1", "2", "3")  
  
  }
  
  
  # ============
  # msm + age
  # ============
  
  if (convergence$msm_age[seed]==2){
    sd <- sqrt(diag(model.msm_age$covmat))[c(4:6, 8:10, 12:14)]
    sd_msm_age <- matrix(sd, nrow = 3, ncol = 3, byrow=T) # are ordered by transition ex 1.cov1,2.cov1,3.cov1,1.cov2,2.cov2..
    
    rownames(sd_msm_age) <- c("cov1", "cov2", "cov3")
    colnames(sd_msm_age) <- c("1", "2", "3")  
    
  } else{
    sd_msm_age <- matrix(NA, nrow = 3, ncol = 3)
    rownames(sd_msm_age) <- c("cov1", "cov2", "cov3")
    colnames(sd_msm_age) <- c("1", "2", "3") 

  }
  
  
  # ========
  # nhm
  # ========
  
  if (convergence$nhm[seed]==2){
    std.err <- diag(solve(model_nhm$hess))^0.5 # ordered like cov1.1 cov1.2 cov1.3
    sd_nhm <- matrix(std.err[7:15], nrow = 3, ncol = 3, byrow=T)
    
    rownames(sd_nhm) <- c("cov1", "cov2", "cov3")
    colnames(sd_nhm) <- c("1", "2", "3")  
    
  } else {
    sd_nhm<- matrix(NA, nrow = 3, ncol = 3)
    rownames(sd_nhm) <- c("cov1", "cov2", "cov3")
    colnames(sd_nhm) <- c("1", "2", "3")
  }
  
  # ============ 
  # imputation 
  # ============ 
  
  if(convergence$imputation[seed]==2){
    sd_imputation <- matrix(0, nrow = 3, ncol = 3)
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
    
    se1 <- sqrt(diag(T1))[3:5]
    se2 <- sqrt(diag(T2))[3:5]
    se3 <- sqrt(diag(T3))[3:5]
    
    sd_imputation <- cbind(se1,se2,se3)
    rownames(sd_imputation) <- c("cov1", "cov2", "cov3")
    colnames(sd_imputation) <- c("1", "2", "3")  
    
  } else {
    sd_imputation <- matrix(NA, nrow = 3, ncol = 3)
    rownames(sd_imputation) <- c("cov1", "cov2", "cov3")
    colnames(sd_imputation) <- c("1", "2", "3")  
  }
  
  sd_tot <- rbind(
    cbind(sd_EO, model = "flexsurv_EO", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(sd_cox, model = "coxph", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(sd_flexsurv, model = "flexsurv", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(sd_msm, model = "msm", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(sd_msm_age, model = "msm_age", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(sd_nhm, model = "nhm", seed = seed, covariate=c("cov1","cov2","cov3")),
    cbind(sd_imputation, model = "imputation", seed = seed, covariate=c("cov1","cov2","cov3"))
  )
  

  return(sd_tot)
}

