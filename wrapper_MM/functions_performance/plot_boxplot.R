plot_boxplot <- function(df_long1,df_long2,df_long3,parameter_name){
  
  df_long1 <- df_long1 %>% filter(parameter == parameter_name)
  df_long2 <- df_long2 %>% filter(parameter == parameter_name)
  df_long3 <- df_long3 %>% filter(parameter == parameter_name)
  
  max_bias <- max(max(df_long1$bias),max(df_long2$bias),max(df_long3$bias))
  min_bias <- 0
  
  plot1 <- ggplot(df_long1, aes(x = model, y = bias, fill = model)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 2) +
    labs(title = paste("Trans 1", parameter_name),
         x = "Model",
         y = "Bias Value") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) + 
    theme_minimal() +
    theme(legend.position = "none",  
          axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
    scale_fill_brewer(palette = "Set3") +
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
    scale_fill_brewer(palette = "Set3") +
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
    scale_fill_brewer(palette = "Set3") +
    ylim(c(min_bias, max_bias)) 
  
  
  plot1 + plot2 + plot3 + plot_layout(ncol = 3)
}