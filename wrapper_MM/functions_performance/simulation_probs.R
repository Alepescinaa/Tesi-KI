covariates_distr <- function(n_obs, meanlog, sdlog, covs){
  cov1 <- rlnorm(n=n_obs, meanlog, sdlog)
  min_val <- min(cov1)
  max_val <- max(cov1)
  cov1 <- (cov1-min_val)/(max_val-min_val)
  cov2 <- rbinom(n=n_obs, size=1, prob=as.numeric(covs[2]))  
  cov3 <- rbinom(n=n_obs, size=1, prob=as.numeric(covs[3]))  
  x1 <- 60
  Intercept <- matrix(1,n_pats,1)
  cov_df <- data.frame(Intercept, x1, cov1, cov2, cov3)
  colnames(cov_df)[1]="(Intercept)"
  return(cov_df=cov_df)   
}

create_input_data<- function(n_pats, meanlog, sdlog, covs){
  res <- covariates_distr(n_pats, meanlog, sdlog, covs )
  patients <- res
  patients$patient_id <- 1:n_pats
  strategies <- data.frame(strategy_id = 1)  
  
  hesim_dat <- hesim_data(strategies = strategies,
                          patients = patients)
  input_data <- hesim::expand(hesim_dat, by = c("strategies", "patients"), times = NULL ) 
  return(input_data[, ])
}

simulation <- function(n_pats, fits, meanlog, sdlog, covs) {
  res <- create_input_data(n_pats, meanlog, sdlog, covs )
  input_data <- res
  entry_matrix <- matrix(0,n_pats,4)
  entry_matrix[,1] <- 1
  colnames(entry_matrix) <- c("x1","cov1","cov2","cov3")
  fic_tmat <- mstate::transMat(x = list(c(2),c(3, 4),c(4),c()), names = c("Zero","Dementia-free","Dementia", "Death"))
  det_model <-  params_surv( coefs = list(est = entry_matrix), dist="fixed")
  det_model$n_samples <- 1
  
  object <- params_surv_list(
    det_model,
    create_params(fits[[1]],n=n_pats, uncertainty="none"),
    create_params(fits[[2]],n=n_pats, uncertainty="none"),
    create_params(fits[[3]],n=n_pats, uncertainty="none")
  )
  
  
  
  dismod <- create_IndivCtstmTrans(
    object,
    input_data,
    trans_mat = fic_tmat,
    clock = "forward",
    start_age = 0
  )
  return(dismod)
}
