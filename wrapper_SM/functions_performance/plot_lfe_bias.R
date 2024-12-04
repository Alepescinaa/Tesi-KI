plot_lfe_bias <- function(dementia){
  if (dementia==0){
    p1 <- ggplot(res_bias_lfe[[1]][[1]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Dementia-free life expectancy (PBS 1 year)",  
           x = "Model",
           y = "Rel Bias") +
      theme_bw() +
      geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")  +
      scale_y_continuous(limits = c(-1,1))
    
    p2 <- ggplot(res_bias_lfe[[2]][[1]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Dementia-free life expectancy (PBS 3 years)",  
           x = "Model",
           y = "Rel Bias") +
      theme_bw() +
      geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")  +
      scale_y_continuous(limits = c(-1,1))
    
    p3 <- ggplot(res_bias_lfe[[3]][[1]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Dementia-free life expectancy (PBS 3-6 years)",  
           x = "Model",
           y = "Rel Bias") +
      theme_bw() +
      geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") +
      scale_y_continuous(limits = c(-1,1))
    
    p4 <- ggplot(res_bias_lfe[[4]][[1]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Dementia-free life expectancy (EHR)",  
           x = "Model",
           y = "Rel Bias") +
      theme_bw() +
      geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")  +
      scale_y_continuous(limits = c(-1,1))
    
    p1 + p2 + p3 + p4 + 
      plot_layout(ncol = 2, nrow = 2, 
                  widths = c(0.6, 0.6), 
                  heights = c(0.8, 0.8)) +
      theme(plot.margin = margin(10, 10, 10, 10))
  } else{
    p1 <- ggplot(mean_estimates_lfe[[1]][[2]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Years of life with dementia (PBS 1 year)",  
           x = "Model",
           y = "Rel Bias") +
      theme_bw() +
      geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")  +
      scale_y_continuous(limits = c(-2,2))
    
    p2 <- ggplot(res_bias_lfe[[2]][[2]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Years of life with dementia (PBS 3 years)",  
           x = "Model",
           y = "Rel Bias") +
      theme_bw() +
      geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")  +
      scale_y_continuous(limits = c(-1,1))
    
    p3 <- ggplot(res_bias_lfe[[3]][[2]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Years of life with dementia (PBS 3-6 years)",  
           x = "Model",
           y = "Rel Bias") +
      theme_bw() +
      geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") +
      scale_y_continuous(limits = c(-1,1))
    
    p4 <- ggplot(res_bias_lfe[[4]][[2]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Years of life with dementia (EHR)",  
           x = "Model",
           y = "Rel Bias") +
      theme_bw() +
      geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")  +
      scale_y_continuous(limits = c(-1,1))
    
    p1 + p2 + p3 + p4 + 
      plot_layout(ncol = 2, nrow = 2, 
                  widths = c(0.6, 0.6), 
                  heights = c(0.8, 0.8)) +
      theme(plot.margin = margin(10, 10, 10, 10))
  }
}