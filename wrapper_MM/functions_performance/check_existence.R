check_existence <- function(n_pats, scheme, seed) {
  
  check <- 0
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
    print(seed_dir)
    files_to_load <- c("model_nhm.RData")
      if (file.exists(files_to_load)) {
        check <- 1
      } else {
        warning(paste("File does not exist:", file))
      }
  } else {
    warning(paste("Seed directory does not exist:", seed_dir))
  }
  return(check)
}
