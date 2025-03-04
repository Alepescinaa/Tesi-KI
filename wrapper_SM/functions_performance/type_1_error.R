type_1_error <- function(significant_covs, data, scheme){
  data <- data[[scheme-1]]
  
  df <- data %>%
    # mutate( beta1=cov1, beta2=cov2, beta3=cov3, cov1=NULL, cov2=NULL, cov3=NULL) %>%
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
    filter(y == 0)  
  
  df_long <- df_long %>%
    mutate(
      power_category = case_when(
        x < 5 ~ "Low (<5%)",
        x < 10 ~ "Moderate (5%-10%)",
        TRUE ~ "High (>10%)"
      ),
    )

  
  return (df_long)
}