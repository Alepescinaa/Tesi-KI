run_imputation <- function(data, data_visits, m, type, cores_imp){
  
  # Create copies of input data for modification
  schemeA_data <- data
  schemeA_visits <- data_visits
  
  # Compute the mid-age between last observed dementia-free age and onset age
  temp <- schemeA_data  
  temp$midage <- (schemeA_data$last_bfo + schemeA_data$onset_age) / 2
  
  # Adjust onset age for individuals with observed dementia onset
  temp$onset_age <- ifelse(temp$onset == 1, temp$midage, temp$onset_age)
  
  # Fit a Cox proportional hazards model for dementia onset
  mod_onset <- coxph(Surv(age, onset_age, onset) ~ cov1 + cov2 + cov3, data = temp)
  onsethaz <- basehaz(mod_onset, center = FALSE)
  onsetpar <- mod_onset$coefficients
  
  # Adjust start time for death model based on onset age
  temp$start <- ifelse(temp$onset == 1, temp$midage, temp$age)
  
  # Fit a Cox proportional hazards model for death
  mod_death <- coxph(Surv(start, death_time, dead) ~ cov1 + cov2 + cov3 + onset, data = temp)
  deathhaz <- basehaz(mod_death, center = FALSE)
  deathpar <- mod_death$coefficients
  
  # Compute hazard increments for onset and death
  tzero <- data.frame(matrix(0, 1, 2))
  names(tzero) <- c("time", "hazard")
  
  onsethaz <- rbind(tzero, onsethaz)
  onsethaz$h <- onsethaz$hazard
  onsethaz$h[nrow(onsethaz)] <- NA
  
  # Compute hazard rate for onset
  for (i in 1:(nrow(onsethaz) - 1)) {
    onsethaz$h[i] <- (onsethaz$hazard[i + 1] - onsethaz$hazard[i]) / (onsethaz$time[i + 1] - onsethaz$time[i])
  }
  
  deathhaz <- rbind(tzero, deathhaz)
  deathhaz$h <- deathhaz$hazard
  deathhaz$h[nrow(deathhaz)] <- NA
  
  # Compute hazard rate for death
  for (i in 1:(nrow(deathhaz) - 1)) {
    deathhaz$h[i] <- (deathhaz$hazard[i + 1] - deathhaz$hazard[i]) / (deathhaz$time[i + 1] - deathhaz$time[i])
  }
  
  remove(tzero)
  
  # Merge cumulative hazards from both models
  combhaz <- merge(onsethaz, deathhaz, by.x = "time", by.y = "time", all.x = TRUE, all.y = TRUE, suffixes = c(".12", ".13"))
  
  # Ensure no missing values in hazard columns
  if (is.na(combhaz$hazard.12[1])) combhaz$hazard.12[1] <- 0
  if (is.na(combhaz$hazard.13[1])) combhaz$hazard.13[1] <- 0
  
  for (i in 2:nrow(combhaz)) {
    if (is.na(combhaz$h.12[i])) {
      combhaz$h.12[i] <- combhaz$h.12[i - 1]
      combhaz$hazard.12[i] <- combhaz$hazard.12[i - 1] + combhaz$h.12[i - 1] * (combhaz$time[i] - combhaz$time[i - 1])
    }
    if (is.na(combhaz$h.13[i])) {
      combhaz$h.13[i] <- combhaz$h.13[i - 1]
      combhaz$hazard.13[i] <- combhaz$hazard.13[i - 1] + combhaz$h.13[i - 1] * (combhaz$time[i] - combhaz$time[i - 1])
    }
  }
  
  # Initialize matrices for dementia status and onset age for imputed datasets
  demstatus <- matrix(nrow = nrow(schemeA_data), ncol = m, rep(schemeA_data$onset, m))
  dementage <- matrix(nrow = nrow(schemeA_data), ncol = m, 0)
  
  # Generate random numbers for imputation
  c1 <- matrix(runif(nrow(schemeA_data) * m), nrow = nrow(schemeA_data), ncol = m)
  c2 <- matrix(runif(nrow(schemeA_data) * m), nrow = nrow(schemeA_data), ncol = m)
  
  # Compute linear predictors for the models
  covar <- as.matrix(schemeA_data[, c("cov1", "cov2", "cov3")])
  lp.12 <- covar %*% as.vector(t(onsetpar[c("cov1", "cov2", "cov3")]))
  lp.13 <- covar %*% as.vector(t(deathpar[c("cov1", "cov2", "cov3")]))
  delta <- deathpar["onset"]
  
  # Iterate over patients to impute dementia onset times
  for (i in 1:nrow(schemeA_data)) {
    patient_data <- schemeA_visits[schemeA_visits$patient_id == i,]
    
    # If dementia is observed, estimate onset time
    if (any(patient_data$onset == 1)) {
      index <- which(patient_data$onset == 1)[1]
      a <- patient_data$age[index - 1]
      b <- patient_data$age[index]
      F <- combhaz[combhaz$time > a & combhaz$time < b,]
      nF <- nrow(F)
      
      if (nF > 0) {
        minP <- F$P[1]
        maxP <- F$P[nF]
        if (minP < maxP) {
          F$q <- (F$P - minP) / (maxP - minP)
          for (j in 1:m) dementage[i, j] <- max(F[F$q <= c1[i, j],]$time)
        } else {
          dementage[i, ] <- (a + b) / 2
        }
      } else {
        dementage[i, ] <- (a + b) / 2
      }
    }
  }
  
  all_fits <- vector(mode = "list", length = m)
  
  plan(multisession, workers = cores_imp)
  
  # Run loop in parallel
  all_fits <- future_lapply(1:m, function(j) {
    temp <- schemeA_data
    
    dem_generated <- which(demstatus[,j] == 1)
    dem_observed <- as.numeric(temp$patient_id[temp$onset == 1])
    dem_to_add <- setdiff(dem_generated, dem_observed)
    
    temp$onset[dem_to_add] <- 1
    temp$onset_age[dem_generated] <- dementage[dem_generated, j]
    
    fit_model(temp, type)  
  })
  
  averaged_params <- averaging_params(all_fits, m)
  
  return(list(averaged_params = averaged_params, all_fits = all_fits))
}
