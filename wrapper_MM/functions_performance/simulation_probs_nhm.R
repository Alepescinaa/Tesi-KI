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

simulation_new <- function(n_pats, fits, meanlog, sdlog, covs, ind) {
  res <- create_input_data(n_pats, meanlog, sdlog, covs )
  input_data <- res
  entry_matrix <- matrix(0,n_pats,4)
  entry_matrix[,1] <- 1
  colnames(entry_matrix) <- c("x1","cov1","cov2","cov3")
  fic_tmat <- mstate::transMat(x = list(c(2),c(3, 4),c(4),c()), names = c("Zero","Dementia-free","Dementia", "Death"))
  det_model <-  params_surv( coefs = list(est = entry_matrix), dist="fixed")
  det_model$n_samples <- 1
  
  if (ind==1){ #nhm
    create_params_mine <- function(fits, n = 1, uncertainty = "none", transition) {
      
      estimates <-  matrix(fits$par, nrow = 3, ncol = 5)
      if (transition %in% 1:3) {
        estimates <- estimates[transition, ]
      }
      shape <- matrix(estimates[2], nrow = 1, ncol = 1, 
                      dimnames = list("est", "(Intercept)"))
      
      
      rate <- matrix(estimates[c(1, 3:5)], nrow = 1, ncol = 4, 
                     dimnames = list("est", c("(Intercept)", "cov1", "cov2", "cov3")))
      
      result <- list(
        coefs = list(
          shape = shape,
          rate = rate
        ),
        dist = "gompertz",
        n_samples = 1
      )
      
      class(result) <- "params_surv"
      return(result)
    }
    
  }
  else if (ind==2){ #msm
    create_params_mine <- function(fits, n = 1, uncertainty = "none", transition) {
      
      estimates <-  matrix(fits$estimates, nrow = 3, ncol = 4)
      if (transition %in% 1:3) {
        estimates <- estimates[transition, ]
      }
      shape <- matrix(0, nrow = 1, ncol = 1, 
                      dimnames = list("est", "(Intercept)"))
      
      
      rate <- matrix(estimates[c(1:4)], nrow = 1, ncol = 4, 
                     dimnames = list("est", c("(Intercept)", "cov1", "cov2", "cov3")))
      
      result <- list(
        coefs = list(
          shape = shape,
          rate = rate
        ),
        dist = "gompertz",
        n_samples = 1
      )
      
      class(result) <- "params_surv"
      return(result)
    }
  }
  else if (ind==3){ #msm+age
    create_params_mine <- function(fits, n = 1, uncertainty = "none", transition) {
      
      estimates <-  matrix(fits$estimates, nrow = 3, ncol = 5)
      if (transition %in% 1:3) {
        estimates <- estimates[transition, ]
      }
      shape <- matrix(estimates[5], nrow = 1, ncol = 1, 
                      dimnames = list("est", "(Intercept)"))
      
      
      rate <- matrix(estimates[c(1:4)], nrow = 1, ncol = 4, 
                     dimnames = list("est", c("(Intercept)", "cov1", "cov2", "cov3")))
      
      result <- list(
        coefs = list(
          shape = shape,
          rate = rate
        ),
        dist = "gompertz",
        n_samples = 1
      )
      
      class(result) <- "params_surv"
      return(result)
    }
  }
  
  
  object <- params_surv_list(
    det_model,
    create_params_mine(fits, transition = 1),
    create_params_mine(fits, transition = 2),
    create_params_mine(fits, transition = 3)
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
