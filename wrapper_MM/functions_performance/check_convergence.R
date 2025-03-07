check_convergence <- function(n_pats, scheme, seed) {
  
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
    files_to_check <- c("model_nhm.RData")
    files_to_load <- c("msm_model.RData", "model_msm_age.RData", "cox_model.RData", "flexsurv_model.RData", "results_imp.RData")
      
    if (!file.exists(files_to_check)) {
      model_nhm <- NULL
      warning(paste("File does not exist:", files_to_check))
      } else {
        load(files_to_check)
      }
    for (i in files_to_load ){
      model <- load(i)
    }
  } else {
    warning(paste("Seed directory does not exist:", seed_dir))
  }
  
  models_imp <- results_imp[[2]]
  
  convergence_results <- tibble( # value to zero no convergence of the algorithm 
    converged_coxph = 1,
    converged_flexsurv = 1,
    converged_nhm = 1,
    converged_msm = 1,
    converged_msm_age = 1,
    converged_imp = 1
  )
  
  if (model_cox$info[[4]] != 0) {
    convergence_results$converged_coxph <- 0
  }
  
  if (any(c(fits_gompertz[[1]]$opt$convergence, 
            fits_gompertz[[2]]$opt$convergence, 
            fits_gompertz[[3]]$opt$convergence) != 0)) {
    convergence_results$converged_flexsurv <- 0
  }
  
  if (model.msm$opt$convergence != 0) {
    convergence_results$converged_msm <- 0
  }

  if (model.msm_age$opt$convergence != 0) {
    convergence_results$converged_msm_age <- 0
  }
  
  if (is.null(model_nhm)) {
    convergence_results$converged_nhm <- 0
  }
  
  for (i in 1:length(models_imp)) {
    if (any(c(models_imp[[i]][[1]]$opt$convergence, 
              models_imp[[i]][[2]]$opt$convergence,
              models_imp[[i]][[3]]$opt$convergence) != 0)) {
      convergence_results$converged_imp <- 0
    }
  }
  
  
  hessian_results <- tibble( # value to zero no definite positive hessian -> no convergence to optimum
    hessian_coxph = 1,
    hessian_flexsurv = 1,
    hessian_nhm = 0,
    hessian_msm = 1,
    hessian_msm_age = 1,
    hessian_imp = 1
  )
  


  if (det(model_cox$var < 1e-7)) { # hessian always computable but we set a threshold to check when it is bad conditioned
    hessian_results$hessian_coxph <- 0
  }
  
  if (any(c(eigen(fits_gompertz[[1]]$opt$hessian)$values, 
            eigen(fits_gompertz[[2]]$opt$hessian)$values, 
            eigen(fits_gompertz[[3]]$opt$hessian)$values) < 0)) {
    hessian_results$hessian_flexsurv <- 0
  }
  
  if (any(eigen(model.msm$opt$hessian)$values < 0)) {
    hessian_results$hessian_msm <- 0
  }
  
  if (any(eigen(model.msm_age$opt$hessian)$values < 0)) {
    hessian_results$hessian_msm_age <- 0
  }
  
  if (!is.null(model_nhm) && model_nhm$singular==FALSE) {
    hessian_results$hessian_nhm <- 1
  }
  
  for (i in 1:length(models_imp)) {
    if (any(c(eigen(models_imp[[i]][[1]]$opt$hessian)$values, 
              eigen(models_imp[[i]][[2]]$opt$hessian)$values, 
              eigen(models_imp[[i]][[3]]$opt$hessian)$values) < 0)) {
      hessian_results$converged_imp <- 0
    }
  }
    
  return(list(convergence_results = convergence_results, hessian_results = hessian_results))
  
}


# Checking model convergence :
# - For nhm: If not reached, it indicates a singular transition probability matrix for some time intervals.
#   nhm estimates transition probabilities by solving differential equations from 0 to t.
#   If P(0,t0) is not invertible, an error occurs. The 'split' option can help but may not always work (LSODA warnings).
#
# - For msm and msm_age: If the covariance matrix of estimated parameters is NULL, confidence intervals won't be computed.
#   This happens when the Hessian is not positive definite, meaning the likelihood’s global maximum was not reached.
#   If optim reports convergence, it met stopping criteria, but parameters may still be improvable.
#   Alternative criteria, initial values, or optimization methods could be explored.
