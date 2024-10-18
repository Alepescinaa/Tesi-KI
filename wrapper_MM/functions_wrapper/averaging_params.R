averaging_params <- function(all_fits,m){
  param_names <- c("shape", "rate", "cov1", "cov2", "cov3" )
  n_params <- length(param_names)
  averaged_params <- matrix(0, nrow = 3, ncol = n_params)
  colnames(averaged_params) <- param_names
  
  for (i in 1:3) {
    for (p in 1:n_params) {
      total_param_value <- 0
      for (j in 1:m) {
        total_param_value <- total_param_value + coef(all_fits[[j]][[i]])[p]
      }
      averaged_params[i, p] <- total_param_value / m
    }
  }
  return(averaged_params)
}