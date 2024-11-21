plot_convergence <- function(scheme, titles){
  df <- combined_cov[[scheme-1]]
  count<- df %>%
    pivot_longer(everything(), names_to = "model", values_to = "convergence") %>%
    group_by(model) %>%
    summarise(
      code_0 = mean(convergence == 0, na.rm = TRUE),
      code_1 = mean(convergence == 1, na.rm = TRUE),
      code_2 = mean(convergence == 2, na.rm = TRUE)
    )%>%
  pivot_longer(cols = starts_with("code_"), names_to = "status", values_to = "percentage")
  
  count <- count %>% filter(status != "code_2")
  count <- count %>%
    mutate(percentage = round(percentage, 3)) 
  
  ggplot(count, aes(x = model, y = percentage, fill = status)) +
    geom_bar(stat = "identity", width = 0.8, position="dodge") + # position = "dodge" se voglio 3 barre affiancate
    scale_fill_manual(values = c(
      "code_0" = "#FF9999",  
      "code_1" = "#99CCFF" 
     # "code_2" = "#99FF99"   
    ), labels = c(
      "code_0" = "Lack of convergence",
      "code_1" = "Suboptimal convergence"
     # "code_2" = "Convergence at Optimum"
    )) +
    labs(title = titles[scheme-1],
         x = "Model",
         y = "Percentage",
         fill = "Convergence status") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "top",  
      legend.title = element_text(size = 9), 
      legend.text = element_text(size = 7),    
      panel.grid.major = element_line(color = "lightgrey", size = 0.3),  
      panel.grid.minor = element_blank()  # No minor grid lines
    )+
    scale_y_continuous(limits = c(0, 0.2)) 
    
  
}