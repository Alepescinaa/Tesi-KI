check_existence <- function(n_pats, scheme, seed) {
  
  nhm_computed <- 0
  converged_msm <- 1
  converged_msm_age <- 1
  
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
    files_to_load <- c("msm_model.RData", "model_msm_age.RData")
      if (file.exists(files_to_check)) {
        nhm_computed <- 1
      } else {
        warning(paste("File does not exist:", files_to_check))
      }
    for (i in files_to_load ){
      model <- load(i)
    }
  } else {
    warning(paste("Seed directory does not exist:", seed_dir))
  }
  
  if (model.msm$opt$convergence!=0)
    converged_msm <- 0
  if (model.msm_age$opt$convergence!=0) 
    converged_msm_age <- 0
  
  return(list(nhm_computed=nhm_computed, converged_msm=converged_msm, converged_msm_age=converged_msm_age ))
}
