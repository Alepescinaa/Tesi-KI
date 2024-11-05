plot_bias <- function (scheme, titles, transition_num){
  
  df <- res_bias[[scheme-1]]
  df_long <- df %>%
    pivot_longer(cols = c("rate", "shape", "cov1", "cov2", "cov3"),
                 names_to = "parameter",
                 values_to = "bias")
  
 
  p <- ggplot(df_long, aes(x = model, y = bias, color = model, shape = transition)) +
    geom_point(size = 3) +
    facet_wrap(~ parameter, scales = "free") +  # Facet by parameter
    labs(title = titles[scheme - 1],  
         x = "Model",
         y = "Bias") +
    theme_bw() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")  
  
  p + geom_line(data = df_long %>% 
                  filter(parameter %in% c("cov1", "cov2", "cov3")), 
                aes(group = transition), 
                size = 0.5, 
                linetype = "dashed", 
                color = "grey", 
                na.rm = TRUE)
  
}
