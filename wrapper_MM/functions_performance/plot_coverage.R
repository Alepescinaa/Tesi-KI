plot_coverage <- function (scheme, titles){
  
  df <- res_cov[[scheme-1]]
  df_long <- df %>%
    pivot_longer(cols = c("rate", "shape", "exp(cov1)", "exp(cov2)", "exp(cov3)"),
                 names_to = "parameter",
                 values_to = "bias")
  df_long <- df_long %>%
    select(-cov1, -cov2, -cov3)
  
  ggplot(df_long, aes(x = model, y = bias, fill = model)) +
    geom_point() +  
    facet_wrap(~ parameter, scales = "free") +
    labs(title = "Bias Comparison Across Models",
         x = "Model",
         y = "Bias") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggplot(df_long, aes(x = model, y = bias, color = model, shape = transition)) +
    geom_point(size = 3) +  # You can adjust the size of the points if needed
    facet_wrap(~ parameter, scales = "free") +
    labs(title = titles[scheme-1],
         x = "Model",
         y = "Bias") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_manual(values = c("coxph" = "red", "flexsurv" = "blue", "imputation" = "green", "msm" = "purple", "msm_age"="black", "nhm"="orange")) +  # Customize colors
    scale_shape_manual(values = c("1"=16, "2"=17, "3"=18)) 
  
}