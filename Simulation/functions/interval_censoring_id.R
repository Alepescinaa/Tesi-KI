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


interval_censoring_id <- function(temp_id, id, scheme, followup, keep_index){
  
  if (scheme == 1) {
    avg_spacing <- 1
    n_visits <- ceiling(followup / avg_spacing)
    obs_scheme <- runif(n_visits - 1, min = avg_spacing - 0.3, max = avg_spacing + 0.3)
  } else if (scheme == 2) {
    avg_spacing <- 3
    n_visits <- ceiling(followup / avg_spacing)
    obs_scheme <- runif(n_visits - 1, min = avg_spacing - 0.5, max = avg_spacing + 0.5)
  } else if (scheme == 3) {
    avg_spacing_low <- 6
    avg_spacing_high <- 3
    n_visits <- ceiling(followup / ((avg_spacing_low + avg_spacing_high) / 2))
    start_age <- temp_id$time_start[temp_id$patient_id == id][1]
    obs_scheme <- observation_scheme(start_age, n_visits = n_visits)
  } else if (scheme == 4) {
    n_visits <- ceiling(runif(1, 1, followup - 1))
    obs_scheme <- runif(n_visits, 1, 8) # All intervals at once
  }
  
  visit_times <- c(0, cumsum(obs_scheme))
  visit_times <- visit_times + temp_id$time_start[1]
  
  num_rows <- nrow(temp_id)
  # identifying the index of the first visit for which the patients is no more observed (death/lost to followup)
  if (num_rows == 1) {
    censor_index <- which(visit_times >= temp_id$death_time)[1]
  } else if (num_rows == 2) {
    temp_id <- temp_id[2,]
    censor_index <- which(visit_times >= temp_id$death_time)[1]
    
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
    temp_id$onset <- ifelse(temp_id$time_start > last_visit, 0, 1)
    onset_index <- which(visit_times >= temp_id$time_start)[1]
    temp_id$onset_age<- ifelse(temp_id$onset == 1, visit_times[onset_index], temp_id$death_time)
    
    # if the patients do not develop dementia in the exactly observed dataset we set onset=0 
  } else {
    temp_id$onset <- 0
    temp_id$onset_age <- temp_id$death_time
  }
  
  # finally we add a column indicating the age at each visit for each patient and update the onset status 
  # when dementia is observed

    n_repeats <- length(visit_times)
    temp_id <- temp_id[rep(1, each = n_repeats), ]
    temp_id$age <- visit_times
    temp_id$visits <- 1:length(visit_times)
  
  
  if (id %in% keep_index){
    temp_id$onset <- ifelse(temp_id$age>=temp_id$onset_age,1,0)
  }
  return(temp_id)
} 
