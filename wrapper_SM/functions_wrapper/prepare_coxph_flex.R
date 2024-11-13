prepare_coxph_flex <- function(scheme_data, n_pats){
  
  original_index <- unique(scheme_data$patient_id)
  scheme_data$last_bfo <- numeric()
  for (i in original_index){
    index <- which(scheme_data$onset[scheme_data$patient_id==i]==1)[1]
    scheme_data$last_bfo[scheme_data$patient_id==i] <- scheme_data$age[scheme_data$patient_id==i][index-1]
  }
  for (i in original_index){
    if(any(scheme_data$onset[scheme_data$patient_id==i]==1))
      scheme_data$onset[scheme_data$patient_id==i] <- 1
  }
  
  scheme_data <- scheme_data[order(scheme_data$patient_id),]
  
  row_id <- scheme_data %>%
    group_by(patient_id) %>%  
    summarise(nrows = n())
  
  scheme_data$patient_id <- rep(1:n_pats, times=as.numeric(row_id$nrows))
  
  scheme_data <- scheme_data %>%
    group_by(patient_id) %>%
    slice(1) %>%
    ungroup()
  
  tmat <- mstate::transMat(x = list(c(2, 3),c(3),c()), names = c("Dementia-free","Dementia", "Death")) 
  n_trans <- max(tmat, na.rm = TRUE)
  
  data_long <- msprep(data = scheme_data, trans = tmat, 
                      time = c(NA, "onset_age", "death_time"), 
                      status = c(NA, "onset", "dead"), 
                      keep = c("age","cov1", "cov2", "cov3"),
                      id="patient_id")
  
  data_long$Tstart[data_long$trans<3] <- data_long$Tstart[data_long$trans<3] +data_long$age[data_long$trans<3]
  data_long$time <- data_long$Tstop-data_long$Tstart
  
  return(data_long)
  
}