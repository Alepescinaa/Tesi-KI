plot_bias_lfe <- function (data, scheme, titles,transition_num){
  
  df_h <- data[[scheme-1]][[1]]
  df_dem <- data[[scheme-1]][[2]]
  
  df_long <- df_h %>%
    pivot_longer(cols = "lfe",
                 names_to = "parameter",
                 values_to = "bias")

  p1 <- ggplot(df_long, aes(x = model, y = bias, color = model)) +
    geom_point(size = 3) +
    facet_wrap(~ parameter, scales = "fixed") +
    labs(title = paste("Mean time in dementia free status ", titles[scheme - 1]),  
         x = "Model",
         y = "Bias") +
    theme_grey() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")  
  
  df_long <- df_dem %>%
    pivot_longer(cols = "lfe",
                 names_to = "parameter",
                 values_to = "bias")
  
  p2 <- ggplot(df_long, aes(x = model, y = bias, color = model)) +
    geom_point(size = 3) +
    facet_wrap(~ parameter, scales = "fixed") +
    labs(title = paste("Mean time in dementia status", titles[scheme - 1]),  
         x = "Model",
         y = "Bias") +
    theme_grey() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")  
  
  p1+p2
  

}