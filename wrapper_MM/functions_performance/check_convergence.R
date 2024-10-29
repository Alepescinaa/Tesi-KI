check_convergence <- function(n_pats, scheme, seed) {
  
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
  
  convergence_results <- tibble( # value to zero no convergence of the algorithm (1 in msm mean reached maxiter)
    converged_coxph = 1,
    converged_flexsurv = 1,
    converged_nhm = 1,
    converged_msm = 1,
    converged_msm_age = 1,
    converged_imp = 1
  )
  
  if (any(c(model_cox[[1]]$info[[4]], model_cox[[2]]$info[[4]], model_cox[[3]]$info[[4]]) != 0)) {
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
  


  if (any(c(det(model_cox[[1]]$var), det(model_cox[[2]]$var), det(model_cox[[3]]$var)) < 1e-4)) {
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
  
  if (!is.null(model_nhm) && all(eigen(model_nhm$hess)$values >0)) {
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


# what are we checking here?
# if the model nhm has been computed
# -> if not it means that the transition probability matrix is singular for some intervals of time
# By default nhm finds the transition probabilities by solving a single initial value problem for each unique covariate pattern. 
# Specifically solves a single system of differential equations starting from 0 to obtain P(0,t) for each t, then uses it to   
# find P(t0,t)=P(0,t0)^-1*P(0,t) for t0>0 
# The method relies on the fact that P(0,t0) is invertible, if the algorithm fails because of the presence of a 
# singular matrix a message error is raised. This problem can be solved by split option in control.
# If a split s is given, P(t0,t) will be find by solving the system of equations for P(t*,t) where t0>t* 
# and t* closest to the s  (that apparently is not always working for me also in those cases I had a warning form 
# LSODA asking for increasing step size )

# if the covariance matrix of the estimated parameter has been computed for msm model and msm_age model
# -> if not we won't compute the confidence intervals. The cov matrix = H^-1 so it is NULL when the Hessian is not definite positive.
# This means that actually the global maximum of the likelihood has not been reached but if optmim reports that convergence is 
# reached it means that the criteria imposed has been satisfied so parameters were not improving anymore. We could try to select
# other criteria or intial values or other method for algorithm...
