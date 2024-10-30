plot_coverage <- function (scheme, titles){
  df <- res_cov[[scheme-1]]
  df_long <- df %>%
    pivot_longer(cols = c("rate", "shape", "cov1", "cov2", "cov3"),
                 names_to = "parameter",
                 values_to = "coverage")
  
  ggplot(df_long, aes(x = model, y = coverage, color = model, shape = transition)) +
    geom_point(size = 3) +  # You can adjust the size of the points if needed
    facet_wrap(~ parameter, scales = "free") +
    labs(title = titles[scheme-1],
         x = "Model",
         y = "Coverage") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_manual(values = c("flexsurv_EO" = "yellow","coxph" = "red", "flexsurv" = "blue", "imputation" = "green", "msm" = "purple", "msm_age"="black", "nhm"="orange")) +  # Customize colors
    scale_shape_manual(values = c("1"=16, "2"=17, "3"=18)) +
    geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", size = 0.5) 
    

}