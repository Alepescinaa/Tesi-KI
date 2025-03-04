onset_count <-  function (data, n_pats){
  temp <- prepare_imputation(data, n_pats)
  m <- 30
  type <- "forward"
  results_imp<- run_imputation_dataset(temp[[1]], temp[[2]], m, type)
  prop_onset <- mean(sapply(1:m, function(i) sum(results_imp[[i]]$onset) / n_pats))
  return(prop_onset)
}
