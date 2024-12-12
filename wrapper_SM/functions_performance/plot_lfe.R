plot_lfe <- function(dementia) {
  if (dementia == 0) {
    p1 <- ggplot(mean_estimates_lfe[[1]][[1]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Dementia-free life expectancy (PBS 1 year)",
           x = "Model",
           y = "Years") +
      theme_bw() +
      geom_hline(yintercept = as.numeric(mean_estimates_lfe[[1]][[1]][1, 2]), 
                 color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 10)) +
      scale_y_continuous(limits = c(15, 23))
    
    p2 <- ggplot(mean_estimates_lfe[[2]][[1]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Dementia-free life expectancy (PBS 3 years)",
           x = "Model",
           y = "Years") +
      theme_bw() +
      geom_hline(yintercept = as.numeric(mean_estimates_lfe[[2]][[1]][1, 2]), 
                 color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 10)) +
      scale_y_continuous(limits = c(15, 23))
    
    p3 <- ggplot(mean_estimates_lfe[[3]][[1]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Dementia-free life expectancy (PBS 3-6 years)",
           x = "Model",
           y = "Years") +
      theme_bw() +
      geom_hline(yintercept = as.numeric(mean_estimates_lfe[[3]][[1]][1, 2]), 
                 color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 10)) +
      scale_y_continuous(limits = c(15, 23))
    
    p4 <- ggplot(mean_estimates_lfe[[4]][[1]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Dementia-free life expectancy (EHR)",
           x = "Model",
           y = "Years") +
      theme_bw() +
      geom_hline(yintercept = as.numeric(mean_estimates_lfe[[4]][[1]][1, 2]), 
                 color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 10)) +
      scale_y_continuous(limits = c(15, 23))
    
    (p1 + p2 + p3 + p4) +
      plot_layout(ncol = 2, nrow = 2, guides = "collect") +
      theme(plot.margin = margin(10, 10, 10, 10),
            legend.position = "right")
    
  } else {
    p1 <- ggplot(mean_estimates_lfe[[1]][[2]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Life expectancy with dementia (PBS 1 year)",
           x = "Model",
           y = "Years") +
      theme_bw() +
      geom_hline(yintercept = as.numeric(mean_estimates_lfe[[1]][[2]][1, 2]), 
                 color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 10)) +
      scale_y_continuous(limits = c(1, 5))
    
    p2 <- ggplot(mean_estimates_lfe[[2]][[2]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Life expectancy with dementia (PBS 3 years)",
           x = "Model",
           y = "Years") +
      theme_bw() +
      geom_hline(yintercept = as.numeric(mean_estimates_lfe[[2]][[2]][1, 2]), 
                 color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 10)) +
      scale_y_continuous(limits = c(1, 4))
    
    p3 <- ggplot(mean_estimates_lfe[[3]][[2]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Life expectancy with dementia (PBS 3-6 years)",
           x = "Model",
           y = "Years") +
      theme_bw() +
      geom_hline(yintercept = as.numeric(mean_estimates_lfe[[3]][[2]][1, 2]), 
                 color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 10)) +
      scale_y_continuous(limits = c(1, 4))
    
    p4 <- ggplot(mean_estimates_lfe[[4]][[2]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Life expectancy with dementia (EHR)",
           x = "Model",
           y = "Years") +
      theme_bw() +
      geom_hline(yintercept = as.numeric(mean_estimates_lfe[[4]][[2]][1, 2]), 
                 color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 10)) +
      scale_y_continuous(limits = c(1, 4))
    
    (p1 + p2 + p3 + p4) +
      plot_layout(ncol = 2, nrow = 2, guides = "collect") +
      theme(plot.margin = margin(10, 10, 10, 10),
            legend.position = "right")
  }
}
