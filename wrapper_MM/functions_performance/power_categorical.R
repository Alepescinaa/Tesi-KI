power_categorical <- function(significant_covs, data, scheme){
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
    filter(y == 1)  
  
  df_long <- df_long %>%
    mutate(
      power_category = case_when(
        x < 50 ~ "Low (<50%)",
        x < 80 ~ "Moderate (50%-80%)",
        TRUE ~ "High (>80%)"
      ),
      covariate = case_when (
        covariate == "cov1"  ~ "beta[1]",
        covariate == "cov2"  ~ "beta[2]",
        covariate == "cov3"  ~ "beta[3]",
        TRUE ~ covariate
      )
    )
  
  df_long1 <- df_long[df_long$transition==1,]
  df_long2 <- df_long[df_long$transition==2,]
  df_long3 <- df_long[df_long$transition==3,]
  
  p1 <- ggplot(df_long1, aes(x = covariate, y = model, fill = power_category)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(
      values = c("Low (<50%)" = "red", "Moderate (50%-80%)" = "yellow", "High (>80%)" = "green"),
      name = "Power Category"
    ) +
    labs(
      title = "Heatmap of Power by Method and Covariate for transition from Dementia-free to Dementia",
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank() 
    )+
    scale_x_discrete(labels = function(x) parse(text = x))
  

  
  p2 <- ggplot(df_long2, aes(x = covariate, y = model, fill = power_category)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(
      values = c("Low (<50%)" = "red", "Moderate (50%-80%)" = "yellow", "High (>80%)" = "green"),
      name = "Power Category"
    ) +
    labs(
      title = "Heatmap of Power by Method and Covariate for transition from Dementia-free to Death",
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.title.x = element_blank(), 
      axis.title.y = element_blank() 
    ) +
    scale_x_discrete(labels = function(x) parse(text = x))
  
  p3 <- ggplot(df_long3, aes(x = covariate, y = model, fill = power_category)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(
      values = c("Low (<50%)" = "red", "Moderate (50%-80%)" = "yellow", "High (>80%)" = "green"),
      name = "Power (%)"
    ) +
    labs(
      title = "Heatmap of Power by Method and Covariate for transition from Dementia to Death",
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.title.x = element_blank(), 
      axis.title.y = element_blank() 
    )+
    scale_x_discrete(labels = function(x) parse(text = x))
  
  return (list(p1,p2,p3))
}