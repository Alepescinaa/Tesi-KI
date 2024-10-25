plot_convergence <- function(scheme, titles){
  df <- convergence_schemes[[scheme-1]]
  percentage_row <- df %>% 
    slice(n()) %>%  
    select(nhm_computed, converged_msm, converged_msm_age, converged_coxph) %>% 
    pivot_longer(everything(), names_to = "model", values_to = "percentage")
  
  ggplot(percentage_row, aes(x = model, y = percentage, fill = model)) +
    geom_bar(stat = "identity") +
    labs(title = titles[scheme-1],
         x = "Model",
         y = "Percentage") +
    theme_minimal() +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) 
  
}