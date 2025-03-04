onset_count <-  function (data, n_pats){
  temp <- prepare_imputation(data, n_pats)
  m <- 30
  type <- "mix"
  error <- F
  tryCatch({
    results_imp<- run_imputation_dataset(temp[[1]], temp[[2]], m, type, cores_nhm)
    
  },
  error = function(e) {
    print(paste("Error during model fitting:", e$message))
    error <<- TRUE
  })
  
  
  if (error) {
    print(paste("No imp convergence"))
    return(NULL)
  }
  
  prop_onset <- mean(sapply(1:m, function(i) sum(results_imp[[i]]$onset) / n_pats))
  return(prop_onset)
}
