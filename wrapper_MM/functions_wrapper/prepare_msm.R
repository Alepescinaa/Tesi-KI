prepare_msm <- function(df){
  df$state <- 1
  df <- df %>%
    group_by(patient_id) %>%
    mutate(
      state = ifelse(onset == 1, 2, state),
      state = ifelse(row_number() == n() & dead == 1, 3, state)
    ) %>%
    ungroup() 
  
  return(df)
}
