plot_bias_rel <- function (scheme, titles,transition_num){
  
  df <- res_bias_rel[[scheme-1]]
  df_long <- df %>%
    pivot_longer(cols = c("rate", "shape", "cov1", "cov2", "cov3"),
                 names_to = "parameter",
                 values_to = "bias") %>%
     filter(transition == transition_num)


  ggplot(df_long, aes(x = model, y = bias, color = model, shape = transition)) +
    geom_point(size = 3) +
    #geom_line(aes(group = transition), size = 0.5, linetype = "dashed", color = "grey", na.rm = TRUE) + 
    facet_wrap(~ parameter, scales = "fixed") +
    labs(title = titles[scheme - 1],  
         x = "Model",
         y = "Relative Bias") +
    theme_grey() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")  
   
 
  
}