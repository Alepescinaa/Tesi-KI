plot_boxplot <- function(df_long1,df_long2,df_long3,parameter_name){
  
  library(RColorBrewer)
  
  df_long1 <- df_long1 %>% filter(parameter == parameter_name)
  df_long2 <- df_long2 %>% filter(parameter == parameter_name)
  df_long3 <- df_long3 %>% filter(parameter == parameter_name)
  
  max_bias <- max(quantile(df_long1$bias, 0.85, na.rm = TRUE),
                  quantile(df_long2$bias, 0.85, na.rm = TRUE),
                  quantile(df_long3$bias, 0.85, na.rm = TRUE))
  min_bias <- 0
  
  color_palette <- brewer.pal(n = 7, name = "Set3")
  model_names <- c("flexsurv_EO", "coxph", "flexsurv", "imputation", "msm", "msm_age", "nhm")
  color_mapping <- setNames(color_palette[1:length(model_names)], model_names)
  
  plot1 <- ggplot(df_long1, aes(x = model, y = bias, fill = model)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 2) +
    labs(title = paste("Trans 1", parameter_name),
         x = "Model",
         y = "Bias Value") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) + 
    theme_minimal() +
    theme(legend.position = "none",  
          axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
    scale_fill_manual(values = color_mapping) +
    ylim(c(min_bias, max_bias)) 
  
  
  plot2 <-ggplot(df_long2, aes(x = model, y = bias, fill = model)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 2) +
    labs(title = paste("Trans 2", parameter_name),
         x = "Model",
         y = "Bias Value") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) + 
    theme_minimal() +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
    scale_fill_manual(values = color_mapping)  +
    ylim(c(min_bias, max_bias)) 
  
  plot3 <- ggplot(df_long3, aes(x = model, y = bias, fill = model)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 2) +
    labs(title = paste("Trans 3", parameter_name),
         x = "Model",
         y = "Bias Value") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) + 
    theme_minimal() +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
    scale_fill_manual(values = color_mapping) +
    ylim(c(min_bias, max_bias)) 
  
  
  plot1 + plot2 + plot3 + plot_layout(ncol = 3)
}