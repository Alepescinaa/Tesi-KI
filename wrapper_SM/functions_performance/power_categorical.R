power_categorical <- function(significant_covs, data, scheme){
  data <- data[[scheme-1]]
  
  df <- data %>%
    mutate(
      model = case_when(
        model == "flexsurv_EO" ~ "0",
        model == "coxph" ~ "a",
        model == "flexsurv" ~ "b",
        model == "imputation" ~ "c",
        TRUE ~ model # Keep other values unchanged if they exist
      )
    )
  
  
  df_merged <- df %>%
    left_join(significant_covs, by = c( "transition"))
  df_merged <- df_merged %>%
    arrange(transition)
  
  df_long <- df_merged %>%
    pivot_longer(
      cols = starts_with("cov"),
      names_to = c("covariate", "type"),
      names_pattern = "(cov[0-9]+)\\.(x|y)"
    ) %>%
    pivot_wider(
      names_from = type,
      values_from = value
    ) %>%
    filter(y == 1)  
  
  df_long <- df_long %>%
    mutate(
      power_category = case_when(
        x < 50 ~ "Low (<50%)",
        x < 80 ~ "Moderate (50%-80%)",
        TRUE ~ "High (>80%)"
      ),
    )

  
  return (df_long)
}