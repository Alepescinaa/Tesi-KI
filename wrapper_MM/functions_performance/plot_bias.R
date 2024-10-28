plot_bias <- function (scheme, titles){
  
  df <- res_bias[[scheme-1]]
  df_long <- df %>%
    pivot_longer(cols = c("rate", "shape", "cov1", "cov2", "cov3"),
                 names_to = "parameter",
                 values_to = "bias")
 
  ggplot(df_long, aes(x = model, y = bias, color = model, shape = transition)) +
    geom_point(size = 3) +
    geom_line(aes(group = transition), size = 0.5, linetype= "dashed", color="grey", na.rm=T) +
    facet_wrap(~ parameter, scales = "free") +  # Facet by parameter
    labs(title = titles[scheme - 1],  # Title from your titles vector
         x = "Model",
         y = "Bias") +
    theme_minimal() +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
    scale_color_manual(values = c("coxph" = "red", 
                                  "flexsurv" = "blue", 
                                  "imputation" = "green", 
                                  "msm" = "purple", 
                                  "msm_age" = "black", 
                                  "nhm" = "orange"))  # Customize colors
    #scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 18))  # Customize shapes
  
}