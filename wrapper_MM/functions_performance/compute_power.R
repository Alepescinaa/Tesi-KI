compute_power <- function(n_pats, scheme, seed, convergence, alpha){
  
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
        print(file)
        print(seed)
        warning(paste("File does not exist:", file))
        model_nhm <- NULL
      }
    }
  } else {
    warning(paste("Seed directory does not exist:", seed_dir))
  }
  
  # ===========
  # EO dataset
  # ===========
  
  power_EO <- matrix(0, nrow = 3, ncol = 3)
  rownames(power_EO) <- c("1->2", "1->3", "2->3")
  colnames(power_EO) <- c("cov1","cov2", "cov3")
  for( i in 1:3){
    z_values <- fits_gompertz_EO[[i]]$res.t[3:5,1]/fits_gompertz_EO[[i]]$res.t[3:5,4]
    p_values <- 2 * (1 - pnorm(abs(z_values)))
    power_EO[i,] <- p_values<=alpha
  }
 


 #  # =========
 #  # coxph
 #  # =========
 
 
  if (convergence$coxph[seed]==2) {
    s <- summary(model_cox)
    p_values <- s$coefficients[,5] #cov1.1,cov2.1,cov3.1... so by row
    power_coxph <- matrix(p_values<=alpha, 3, 3, byrow = T)
    rownames(power_coxph) <- c("1->2", "1->3", "2->3")
    colnames(power_coxph) <- c("cov1","cov2", "cov3")
    power_coxph <- ifelse(power_coxph, 1, 0)
  }else{
    power_coxph <- matrix(NA, nrow = 3, ncol = 3)
    rownames(power_coxph) <- c("1->2", "1->3", "2->3")
    colnames(power_coxph) <- c("cov1","cov2", "cov3")
    
    }
 #  
 #  # ============
 #  # flexsurv
 #  # ============

  if (convergence$flexsurv[seed]==2){
    power_flex <- matrix(0, nrow = 3, ncol = 3)
    rownames(power_flex) <- c("1->2", "1->3", "2->3")
    colnames(power_flex) <- c("cov1","cov2", "cov3")
    for( i in 1:3){
      z_values <- fits_gompertz[[i]]$res.t[3:5,1]/fits_gompertz[[i]]$res.t[3:5,4]
      p_values <- 2 * (1 - pnorm(abs(z_values)))
      power_flex[i,] <- p_values<=alpha
    }
  }else{
    power_flex <- matrix(NA, nrow = 3, ncol = 3)
    rownames(power_flex) <- c("1->2", "1->3", "2->3")
    colnames(power_flex) <- c("cov1","cov2", "cov3")
  }
  gc()

 #  # ============
 #  # msm
 #  # ============

  if (convergence$msm[seed]==2){
    z_values <- model.msm$estimates[4:12]/sqrt(diag(model.msm$covmat)[4:12]) #cov1.1,cov1.2,cov1.3... so by col
    p_values <- 2 * (1 - pnorm(abs(z_values)))
    power_msm <- matrix(p_values<=alpha, 3, 3, byrow = F)
    rownames(power_msm) <- c("1->2", "1->3", "2->3")
    colnames(power_msm) <- c("cov1","cov2", "cov3")
    power_msm <- ifelse(power_msm, 1, 0)
    
  } else {
    power_msm <- matrix(NA, nrow = 3, ncol = 3)
    rownames(power_msm) <- c("1->2", "1->3", "2->3")
    colnames(power_msm) <- c("cov1","cov2", "cov3")
 }

 #  # ============
 #  # msm + age
 #  # ============

  if (convergence$msm_age[seed]==2){
    z_values <- model.msm_age$estimates[4:12]/sqrt(diag(model.msm_age$covmat)[4:12]) #cov1.1,cov1.2,cov1.3... so by col
    p_values <- 2 * (1 - pnorm(abs(z_values)))
    power_msm_age <- matrix(p_values<=alpha, 3, 3, byrow = F)
    rownames(power_msm_age) <- c("1->2", "1->3", "2->3")
    colnames(power_msm_age) <- c("cov1","cov2", "cov3")
    power_msm_age <- ifelse(power_msm_age, 1, 0)
  } else {
    power_msm_age <- matrix(NA, nrow = 3, ncol = 3)
    rownames(power_msm_age) <- c("1->2", "1->3", "2->3")
    colnames(power_msm_age) <- c("cov1","cov2", "cov3")
  }

  gc()
  
 #  # ============
 #  # nhm
 #  # ============

  if (convergence$nhm[seed]==2){
    est <- model_nhm$par[7:15]
    cov_matrix <- diag(solve(model_nhm$hess))[7:15]
    z_values <- est/sqrt(cov_matrix)
    p_values <- 2 * (1 - pnorm(abs(z_values)))
    power_nhm<- matrix(p_values<=alpha, 3, 3, byrow = F)  #cov1.1,cov1.2,cov1.3... so by col
    power_nhm <- ifelse(power_nhm, 1, 0)
    rownames(power_nhm) <- c("1->2", "1->3", "2->3")
    colnames(power_nhm) <- c("cov1","cov2", "cov3")
    
  } else{
    power_nhm<- matrix(NA, 3, 3)  
    rownames(power_nhm) <- c("1->2", "1->3", "2->3")
    colnames(power_nhm) <- c("cov1","cov2", "cov3")
  }

  gc()
  #
 #  # ============
 #  # imputation
 #  # ============

  ## correggere p-value 
  if (convergence$imputation[seed]==2){
    power_imp<- matrix(0, 3, 3)  
    rownames(power_imp) <- c("1->2", "1->3", "2->3")
    colnames(power_imp) <- c("cov1","cov2", "cov3")
    
    m <- length(results_imp[[2]])
    all_fits <- results_imp[[2]]
    avg_parameters <- results_imp[[1]]
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
    
    cov_matrices <- list(diag(T1[3:5,3:5]), diag(T2[3:5,3:5]), diag(T3[3:5,3:5]))
    est <- list(avg_parameters[1,3:5],avg_parameters[2,3:5],avg_parameters[3,3:5])
  
    for( i in 1:3){
      z_values <- est[[i]]/sqrt(cov_matrices[[i]])
      p_values <- 2 * (1 - pnorm(abs(z_values)))
      power_imp[i,] <- p_values<=alpha
    }
    } else{
    power_imp<- matrix(NA, 3, 3)  
    rownames(power_imp) <- c("1->2", "1->3", "2->3")
    colnames(power_imp) <- c("cov1","cov2", "cov3")

  }


  power_tot <- rbind(
    cbind(power_EO, model = "flexsurv_EO", seed=seed, transition = c(1,2,3)),
    cbind(power_coxph, model = "coxph", seed = seed , transition = c(1,2,3)),
    cbind(power_flex, model = "flexsurv", seed = seed, transition = c(1,2,3)),
    cbind(power_msm, model = "msm", seed = seed, transition = c(1,2,3)),
    cbind(power_msm_age, model = "msm_age", seed = seed,transition = c(1,2,3)),
    cbind(power_nhm, model = "nhm", seed = seed, transition = c(1,2,3)),
    cbind(power_imp, model = "imputation", seed = seed, transition = c(1,2,3))
  )

 
  return (power_tot)
}

