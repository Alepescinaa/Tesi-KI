prepare_SmoootHazard <- function(scheme_data, n_pats){
  
  original_index <- unique(scheme_data$patient_id)
  scheme_data$last_bfo <- numeric()
  for (i in original_index){
    index <- which(scheme_data$onset[scheme_data$patient_id==i]==1)[1]
    scheme_data$last_bfo[scheme_data$patient_id==i] <- scheme_data$age[scheme_data$patient_id==i][index-1]
    if (all(is.na(scheme_data$last_bfo[scheme_data$patient_id==i]))){
      last_visit <- max(scheme_data$visits[scheme_data$patient_id==i])
      patient_data <- scheme_data[scheme_data$patient_id == i, ]
      age_at_last_visit <- patient_data$age[patient_data$visits == last_visit]
      scheme_data$last_bfo[scheme_data$patient_id==i] <- age_at_last_visit
      scheme_data$onset_age[scheme_data$patient_id==i] <- age_at_last_visit
    }
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
  
  scheme_data <- scheme_data %>%
    mutate(l=last_bfo, r=onset_age, last_bfo=NULL, onset_age=NULL)
  
 
  return(scheme_data)
  
}