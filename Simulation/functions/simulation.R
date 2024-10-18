covariates_distr <- function(n_obs, meanlog, sdlog){
  cov1 <- rlnorm(n=n_obs, meanlog, sdlog)
  min_val <- min(cov1)
  max_val <- max(cov1)
  cov1 <- (cov1-min_val)/(max_val-min_val)
  cov2 <- rbinom(n=n_obs, size=1, prob=0.15)   
  cov3 <- rbinom(n=n_obs, size=1, prob=0.30)  
  x1 <- runif(n_obs, min = 60, max = 96)
  Intercept <- matrix(1,n_pats,1)
  cov_df <- data.frame(Intercept, x1, cov1, cov2, cov3)
  colnames(cov_df)[1]="(Intercept)"
  return(list(cov_df=cov_df, extremes=c(min_val,max_val)))    
}

create_input_data<- function(n_pats, meanlog, sdlog){
  res <- covariates_distr(n_pats, meanlog, sdlog )
  patients <- res[[1]]
  extremes <- res[[2]]
  patients$patient_id <- 1:n_pats
  strategies <- data.frame(strategy_id = 1)  
  
  hesim_dat <- hesim_data(strategies = strategies,
                          patients = patients)
  input_data <- hesim::expand(hesim_dat, by = c("strategies", "patients"), times = NULL ) 
  return(list(input_data[, ], extremes))
}

simulation <- function(n_pats, fits, model_type, meanlog, sdlog) {
  res <- create_input_data(n_pats, meanlog, sdlog )
  input_data <- res[[1]]
  extremes <- res[[2]]
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
    clock = model_type,
    reset_states = 3,
    start_age = 0
  )
  return(list(dismod, extremes))
}
