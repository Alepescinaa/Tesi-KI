level_convergence <- function(scheme){
  combined_cov <- data.frame()
  for (i in 1:4) {
    for (j in 1:100){
      if(convergence_schemes[[scheme-1]][j,i]==0)
        combined_cov[j,i] <- 0
      else if (convergence_schemes[[scheme-1]][j,i]==1 & hessian_schemes[[scheme-1]][j,i]==0 )
        combined_cov[j,i] <- 1
      else if (convergence_schemes[[scheme-1]][j,i]==1 & hessian_schemes[[scheme-1]][j,i]==1 )
        combined_cov[j,i] <- 2
    }
  }
  colnames(combined_cov) <- c("coxph", "flexsurv", "nhm", "imputation")
  return(combined_cov)
  
}