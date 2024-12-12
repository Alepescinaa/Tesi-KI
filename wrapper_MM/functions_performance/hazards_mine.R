hazards_mine <- function(x, b.covariates, no.years, trans = NULL, CI = FALSE, 
         age.shift = 0) {
  
  model <- x  # Fitted MSM model
  
  # Handle transition matrix (convert vector to matrix if needed)
  if (is.vector(trans)) {
    trans <- matrix(trans, 1, 2)
  }
  
  # Create the transition matrix if not provided
  if (!is.null(trans)) {
    ntrans <- nrow(trans)  # Number of transitions
  } else {
    Q <- model$qmodel$imatrix  # Transition rate matrix
    ntrans <- sum(Q)  # Count the number of transitions
    trans <- matrix(NA, ntrans, 2)  # Initialize matrix for transitions
    index <- 1
    for (i in 1:nrow(Q)) {
      for (j in 1:nrow(Q)) {
        if (Q[i, j]) {
          trans[index, ] <- c(i, j)  # Add transitions to matrix
          index <- index + 1
        }
      }
    }
  }
  
  # Number of states in the model
  nstates <- nrow(model$Qmatrices$baseline)
  
  # Create a null transition matrix with zeros on the diagonal
  Q.null <- matrix(as.numeric(model$Qmatrices$baseline != 0), nstates, nstates)
  diag(Q.null) <- 0
  ntrans <- sum(Q.null)  # Count the number of transitions
  
  # Number of covariates (including age and other covariates)
  ncovs <- 1 + length(b.covariates)
  nbeta <- max(which(names(model$estimates) == "qcov"))
  
  if (nbeta != (ntrans * ncovs)) {
    stop("\nAll covariates in msm model should be specified.\n\n")
  }
  
  # Extract the age variable from the covariates
  b.age <- b.covariates$age
  
  # Create a sequence of ages for computation
  age.grid <- seq(b.age, b.age + no.years, by = 1/12)
  
  # Initialize an empty list to store the hazards
  hazards <- list()
  
  # Loop over each transition to calculate the hazards
  for (i in 1:ntrans) {
    haz <- rep(NA, length(age.grid))
    hazLB <- rep(NA, length(age.grid))
    hazUB <- rep(NA, length(age.grid))
    
    # For each age, calculate the hazard for the current transition
    for (j in 1:length(age.grid)) {
      if (ncovs == 2) {
        covariates <- list(age = age.grid[j])
      } else {
        covariates <- c(age = age.grid[j], b.covariates[2:length(b.covariates)])
      }
      
      # Compute the hazard for the current transition (without confidence intervals)
      haz[j] <- qmatrix.msm(model, covariates = covariates, ci = "none")[trans[i, 1], trans[i, 2]]
      
      if (CI) {
        # Compute confidence intervals if requested
        Q.CI <- qmatrix.msm(model, covariates = covariates, ci = "delta")
        hazLB[j] <- Q.CI$L[trans[i, 1], trans[i, 2]]
        hazUB[j] <- Q.CI$U[trans[i, 1], trans[i, 2]]
      }
    }
    
    # Store the hazards (and optionally the confidence intervals)
    hazard_data <- list(hazard = haz)
    if (CI) {
      hazard_data$hazardLB <- hazLB
      hazard_data$hazardUB <- hazUB
    }
    
    # Name the transition for easier identification later
    hazard_name <- paste("Transition", trans[i, 1], "to", trans[i, 2])
    hazards[[hazard_name]] <- hazard_data
  }
  
  # Return the computed hazards (and confidence intervals if requested)
  return(list(hazards, age.grid))
}
