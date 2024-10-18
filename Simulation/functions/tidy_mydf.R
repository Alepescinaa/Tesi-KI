tidy_mydf <- function(df){
  # df$to[df$from==2 & df$to==2] <- 3
  df <- df %>%
    group_by(patient_id) %>%
    slice(ifelse(n() > 1, 2, 1)) %>%
    ungroup()
  df <- df %>%
    mutate(
      death =  1 - censored, 
      onset = ifelse(is.na(onset), 0, onset),        
      onset_age = ifelse(is.na(onset_age), death_time, onset_age)  
    )
  df <- df[,c(1,9:12)]
  return(df)
}