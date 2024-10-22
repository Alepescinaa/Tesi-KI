fit_model <- function(temp, type){
  
  
  tmat <- mstate::transMat(x = list(c(2, 3), c(3), c()), names = c("Dementia-free", "Dementia", "Death"))
  
  data_long <- msprep(data = temp, trans = tmat,
                      time = c(NA, "onset_age", "death_time"),
                      status = c(NA, "onset", "dead"),
                      keep = c("age", "cov1", "cov2", "cov3"),
                      id = "patient_id")
  
 
  data_long$Tstart[data_long$trans < 3] <- data_long$Tstart[data_long$trans < 3] + data_long$age[data_long$trans < 3]
  data_long$time <- data_long$Tstop - data_long$Tstart
  
  n_trans <- max(tmat, na.rm = TRUE)
  fits_wei <- vector(mode = "list", length = n_trans)
  

  
  if (type == "forward") {
    for (i in 1:3) {
      fits_wei[[i]] <- flexsurvreg(Surv(Tstart, Tstop, status) ~ cov1 + cov2 + cov3, 
                                   data = subset(data_long, trans == i), 
                                   dist = "gompertz")
    }
  } else if (type == "mix") { 
    for (i in 1:2) {
      fits_wei[[i]] <- flexsurvreg(Surv(Tstart, Tstop, status) ~ cov1 + cov2 + cov3, 
                                   data = subset(data_long, trans == i), 
                                   dist = "gompertz")
    }
    fits_wei[[3]] <- flexsurvreg(Surv(time, status) ~ cov1 + cov2 + cov3, 
                                 data = subset(data_long, trans == 3), 
                                 dist = "gompertz")
  }
  
  

  return(fits_wei)
}