plot_se <- function (data,scheme, titles){
  df <- data[[scheme-1]]
  df <- as.data.frame(df)
  desired_order <- c("flexsurv_EO", "coxph", "flexsurv", "msm", "msm_age", "nhm", "imputation")
  df$model <- factor(df$model, levels = desired_order, ordered = TRUE)
  df <- df[order(df$model), ]
  df_long <- df %>%
    pivot_longer(cols = c("1", "2", "3"),
                 names_to = "transition",
                 values_to = "se_value") 

  ggplot(df_long, aes(x = model, y = se_value, color = model, shape = transition)) +
    geom_point(size=3) +
    facet_grid(transition~ covariate, scales = "fixed") +  
    coord_cartesian(ylim = c(0, 3))+
    labs(title = titles[scheme - 1],  
         x = "Model",
         y = "Mean SE") +
    theme_bw() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  
 
  
}