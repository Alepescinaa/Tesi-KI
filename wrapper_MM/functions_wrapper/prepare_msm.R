prepare_msm <- function(df){
  df$state <- 1
  df <- df %>%
    group_by(patient_id) %>%
      do({
        if (nrow(.) == 1) {
          duplicated_row <- bind_rows(., .)
          duplicated_row$age[2] <- duplicated_row$death_time[1]
          duplicated_row$visits[2] <- 2
          duplicated_row
        } else {
          .
        }
      }) %>%
      ungroup()
  df <- df %>%
    group_by(patient_id) %>%
    mutate(
      state = ifelse(onset == 1, 2, state),
      state = ifelse(row_number() == n() & dead == 1, 3, state),
      state = ifelse(row_number() == n() & dead == 0 & onset == 0, 99, state)
    ) %>%
    ungroup() 
  
  return(df)
}
