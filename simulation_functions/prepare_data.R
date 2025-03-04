prepare_data <- function(dataset){
  dataset$last_obs <-pmax(dataset$age_wave1, 
                          dataset$age_wave2, 
                          dataset$age_wave3, 
                          dataset$age_wave4, 
                          dataset$age_wave5, 
                          dataset$age_wave6,
                          na.rm = TRUE) 
  
  
  end_of_study <- max(dataset$Date_wave6, na.rm=T)
  dataset$Date_wave6[dataset$Date_wave6==end_of_study] <- end_of_study+1
  dataset$age_end_of_study<- dataset$age_wave6 + interval(dataset$Date_wave6, end_of_study) / years(1)
  
  
  dataset$censoring <- ifelse(dataset$last_obs<78, dataset$last_obs + 6, dataset$last_obs + 3) 
  dataset$censoring <- ifelse(!is.na(dataset$dementia_wave6), dataset$age_end_of_study, dataset$censoring)
  
  dataset$death_time <- ifelse(dataset$dead == 1, dataset$age_death, dataset$censoring)
  
  age_columns <- c("age_wave1", "age_wave2", "age_wave3", "age_wave4", "age_wave5", "age_wave6")
  dementia_columns <- c("dementia_wave1", "dementia_wave2", "dementia_wave3", 
                        "dementia_wave4", "dementia_wave5", "dementia_wave6")
  
  dataset <- dataset %>%
    mutate(across(dementia_columns, ~ as.factor(as.character(.))))
  
  dataset$onset_age <- apply(dataset, 1, function(row) {
    onset_index <- which(row[dementia_columns] == 2)[1] #in data raw dementia presence identified by 2
    
    if (!is.na(onset_index)) { 
      return(row[age_columns[onset_index]]) 
    } else {
      return(NA) 
    }
  })
  
  dataset$onset <- rep(0,nrow(dataset))
  dataset <- dataset %>%
    mutate(onset = if_else(!is.na(onset_age), 1, 0))
  
  dataset <- dataset %>%
    mutate(onset_age = if_else(is.na(onset_age), death_time, as.numeric(onset_age))) 
  return(dataset)
}


#mstate needs as input 
# - onset : binary variable indicating who developed dementia
# - dead : binary variable indicating whose death was observed
# - onset_age : time for which transition to dementia is at risk
# - death_time : time for which transition to death is at risk

#censoring patients death-free
# - those lost to follow up censored at the time of the next expected visit
# - those death-free upt to wave 6 censored at the end of the study (censoring patient observed exactly at the end of study one day later)

#defining death_time
# -age_death if dead=1
# -censoring value if dead=0

#defining onset_age
# -first visit at which dementia was observed if dementia=1
# -death_time if dementia=0
