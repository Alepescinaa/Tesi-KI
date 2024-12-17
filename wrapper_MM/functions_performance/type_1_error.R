type_1_error <- function(significant_covs, data, scheme){
  data <- data[[scheme-1]]
  
  df <- data %>%
    # mutate( beta1=cov1, beta2=cov2, beta3=cov3, cov1=NULL, cov2=NULL, cov3=NULL) %>%
    mutate(
      model = case_when(
        model == "flexsurv_EO" ~ "0",
        model == "coxph" ~ "a",
        model == "flexsurv" ~ "b",
        model == "msm" ~ "c",
        model == "msm_age" ~ "d",
        model == "nhm" ~ "e",
        model == "imputation" ~ "f",
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
      covariate = case_when (
        covariate == "cov1"  ~ "beta1",
        covariate == "cov2"  ~ "beta2",
        covariate == "cov3"  ~ "beta3",
        TRUE ~ covariate
      )
    )
  
  df_long1 <- df_long[df_long$transition==1,]
 
  df_long3 <- df_long[df_long$transition==3,]
  
  p1 <- ggplot(df_long1, aes(x = covariate, y = model, fill = power_category)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(
      values = c("Low (<5%)" = "green", "Moderate (5%-10%)" = "yellow", "High (>10%)" = "red"),
      name = "Type I Error Category"
    ) +
    labs(
      x = "Transition Dementia-free -> Dementia",
      y = "Model"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, hjust = 0.5)
    )
  

  
  p3 <- ggplot(df_long3, aes(x = covariate, y = model, fill = power_category)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(
      values = c("Low (<5%)" = "green", "Moderate (5%-10%)" = "yellow", "High (>10%)" = "red"),
      name = "Type I Error Category"
    ) +
    labs(
      x = "Transition Dementia -> Death",
      y = "Model"

    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, hjust = 0.5)
    )
  
  combined_plot <- (p1 + p3) +
    plot_layout(guides = "auto") + 
    plot_annotation(
      title = "Heatmap of Type I Error by Method and Covariate",
      theme = theme(
        plot.title = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(size = 8, hjust = 0.5)
      )
    )
  
  return (combined_plot)
}