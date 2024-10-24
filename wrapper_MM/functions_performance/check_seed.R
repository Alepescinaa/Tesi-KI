check_seed <- function(scheme_dir, seed) {
  check <- 0
  seed_dir <- file.path(scheme_dir, paste0("seed_", seed))
  
  if (dir.exists(seed_dir)) {
    model_file <- file.path(seed_dir, "model_nhm.RData")
    
    if (!file.exists(model_file)) {
      warning(paste("model_nhm does not exist for seed:", seed))
    } else {
      check <- 1
    }
  } else {
    warning(paste("Seed directory does not exist:", seed_dir))
  }
  
  return(check)
}