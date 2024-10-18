prepare_imputation <- function(data){
  scheme_visits <- data
  scheme_data <- data
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
  
  row_id <- scheme_visits %>%
    group_by(patient_id) %>%  
    summarise(nrows = n())
  
  scheme_visits$patient_id <- rep(1:n_pats, times=as.numeric(row_id$nrows))
  scheme_data$patient_id <- rep(1:n_pats, times=as.numeric(row_id$nrows))
  
  scheme_data <- scheme_data %>%
    group_by(patient_id) %>%
    slice(1) %>%
    ungroup()
  
  return(list(scheme_data=scheme_data, scheme_visits=scheme_visits ))
  
}