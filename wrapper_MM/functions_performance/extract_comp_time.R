extract_comp_time <- function(n_pats, scheme, seed){
  
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
    
    files_to_load <- c("computational_time.RData")
    load(files_to_load)
  } else {
    warning(paste("Seed directory does not exist:", seed_dir))
  }
  
  return(comp_time)
}
