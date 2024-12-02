plot_width <- function (data,scheme, titles){
  df <- data[[scheme-1]]
  df <- as.data.frame(df)
  desired_order <- c("flexsurv_EO", "coxph", "flexsurv", "imputation")
  df$model <- factor(df$model, levels = desired_order, ordered = TRUE)
  df <- df[order(df$model), ]
  df_long <- df %>%
    pivot_longer(cols = c("1", "2", "3"),
                 names_to = "transition",
                 values_to = "ic_value") 

  ggplot(df_long, aes(x = model, y = ic_value, color = model, shape = transition)) +
    geom_line(size=3) +
    #geom_point(size=3)+
    facet_grid(transition~ covariate, scales = "free") +  
    labs(title = titles[scheme - 1],  
         x = "Model",
         y = "Mean width") +
    theme_bw() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  
 
  
}