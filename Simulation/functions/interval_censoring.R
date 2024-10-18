# function for generating observation scheme of type 3, i.e every 6 year for patients with age < 78
# and every 3 years for patients with age > 78

observation_scheme <- function(start_age, n_visits = 6, age_threshold = 78, avg_spacing_low = 6, avg_spacing_high = 3) {
  
  # Create logical vector to check if current age is below or above threshold
  is_below_threshold <- (start_age + cumsum(rep(avg_spacing_low, n_visits - 1)) - avg_spacing_low) < age_threshold
  
  # Use vectorized runif to generate intervals based on the threshold condition
  obs_scheme <- ifelse(is_below_threshold,
                       runif(n_visits - 1, min = avg_spacing_low - 0.5, max = avg_spacing_low + 0.5),
                       runif(n_visits - 1, min = avg_spacing_high - 0.5, max = avg_spacing_high + 0.5))
  
  return(obs_scheme)
}


interval_censoring <- function(temp, id, scheme, followup, keep_index){
  
  
  if (scheme==1){ # every year
    avg_spacing <- 1
    n_visits <- ceiling(followup/avg_spacing)
    obs_scheme <- runif(n_visits - 1, min = avg_spacing - 0.3, max = avg_spacing + 0.3)
  }
  if (scheme==2){ # every 3 years
    avg_spacing <- 3
    n_visits <- ceiling(followup/avg_spacing)
    obs_scheme <- runif(n_visits - 1, min = avg_spacing - 0.5, max = avg_spacing + 0.5)
  }
  if (scheme==3) # snac-k scheme
  {
    avg_spacing_low <- 6
    avg_spacing_high <- 3
    n_visits <- ceiling(followup/((avg_spacing_low + avg_spacing_high)/2))
    start_age <- temp$time_start[temp$patient_id==id][1]
    obs_scheme <- observation_scheme(start_age, n_visits = n_visits)
  }
  if (scheme==4) # Uk Biobank scheme 
  {
    n_visits <- ceiling(runif(1,1,followup-1)) #set random number of visits for each patients
    obs_scheme <- numeric()
    for (i in 1:n_visits)
      obs_scheme[i]<- runif(1,1,8) #set interval between two visits as random number extracted each time between (1,8) years
  }
  visit_times <- c(0, cumsum(obs_scheme))
  visit_times <- visit_times + temp$time_start[temp$patient_id==id][1]
  num_rows <- sum(temp$patient_id == id)
  # identifying the index of the first visit for which the patients is no more observed (death/lost to followup)
  if (num_rows == 1) {
    censor_index <- which(visit_times >= temp$death_time[temp$patient_id == id])[1]
  } else if (num_rows == 2) {
    censor_index <- which(visit_times >= temp$death_time[temp$patient_id == id][2])[1]
  }
  # if no index found it means that the patient is observed at all waves
  # we identify the age at last observed visit
  if (!is.na(censor_index)) {
    censor_index <- ifelse(censor_index > 1, censor_index, 2)
    visit_times <- visit_times[1:(censor_index - 1)]
    last_visit <- visit_times[censor_index - 1]
  } else {
    last_visit <- visit_times[n_visits]
  }
  # if the patients develops dementia in the exactly observed dataset we set onset=1 in the censored dataset
  # only if the (real) dementia onset is observed before last attended visit, onset=0 otherwise
  # in the former case we set onset_age at the first visit for which age > (real) dementia onset age in the exactly observed dataset
  # in the latter we set it at censoring time (death/lost to follow up)
  if (id %in% keep_index) {
    temp$onset[temp$patient_id == id][2] <- ifelse(temp$time_start[temp$patient_id == id][2] > last_visit, 0, 1)
    onset_index <- which(visit_times >= temp$time_start[temp$patient_id == id][2])[1]
    temp$onset_age[temp$patient_id == id][2] <- ifelse(temp$onset[temp$patient_id == id][2] == 1, 
                                                       visit_times[onset_index], 
                                                       temp$death_time[temp$patient_id == id][2])
    
    # if the patients do not develop dementia in the exactly observed dataset we set onset=0 
  } else {
    temp$onset[temp$patient_id == id] <- 0
    temp$onset_age[temp$patient_id == id] <- temp$death_time[temp$patient_id == id]
  }
  
  # finally we add a column indicating the age at each visit for each patient and update the onset status 
  # when dementia is observed
  if (!(id %in% keep_index)){
    patient_data <- temp[temp$patient_id == id, ]
    n_repeats <- length(visit_times)
    repeated_patient_data <- patient_data[rep(1, each = n_repeats), ]
    temp <- temp[temp$patient_id != id, ]  
    temp <- rbind(temp, repeated_patient_data)
    temp$age[temp$patient_id == id] <- visit_times
    temp$visits[temp$patient_id == id] <- 1:length(visit_times)
  }
  
  
  if (id %in% keep_index){
    patient_data <- temp[temp$patient_id == id, ][2]
    n_repeats <- length(visit_times)
    repeated_patient_data <- patient_data[rep(1, each = n_repeats), ]
    temp <- temp[temp$patient_id != id, ]  
    temp <- rbind(temp, repeated_patient_data)
    temp$age[temp$patient_id == id] <- visit_times
    temp$visits[temp$patient_id == id] <- 1:length(visit_times)
    temp$onset <- ifelse(temp$age>=temp$onset_age,1,0)
  }
  return(temp)
}