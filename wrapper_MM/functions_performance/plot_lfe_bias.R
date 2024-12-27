
plot_lfe_bias <- function(dementia) {
  if (dementia == 0) {
    
    create_plot <- function(data, title, y_label = "Rel Bias", y_limits = c(-1, 1)) {
      ggplot(data, aes(x = model, y = lfe, color = model)) +
        geom_point(size = 4) +
        labs(
          title = title,
          x = "Model",
          y = y_label
        ) +
        theme_bw() +
        geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 9)
        ) +
        scale_y_continuous(limits = y_limits)
    }
    ggplot(data, aes(x = model, y = lfe_mean, fill = as.factor(dem))) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
      geom_errorbar(aes(ymin = lfe_lower_ci, ymax = lfe_upper_ci), 
                    position = position_dodge(width = 0.7), width = 0.25) +
      labs(
        title = "Relative Bias of Time Spent in Status",
        x = "Model",
        y = "Relative Bias",
        fill = "State (dem)"
      ) +
      scale_y_continuous(limits = c(-1.2, 0)) +
      scale_fill_manual(values = c("0" = "skyblue", "1" = "orange"),
                        labels = c("Without Dementia", "With Dementia")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5)
      )
    (p1 + p2 + p3 + p4) +
      plot_layout(ncol = 2, nrow = 2, guides = "collect") +
      theme(plot.margin = margin(10, 10, 10, 10),
            legend.position = "right")
    
  } else {
    p1 <- ggplot(res_bias_lfe[[1]][[2]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Years of life with dementia (PBS 1 year)",
           x = "Model",
           y = "Rel Bias") +
      theme_bw() +
      geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 9)) +
      scale_y_continuous(limits = c(-2, 2))
    
    p2 <- ggplot(res_bias_lfe[[2]][[2]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Years of life with dementia (PBS 3 years)",
           x = "Model",
           y = "Rel Bias") +
      theme_bw() +
      geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 9)) +
      scale_y_continuous(limits = c(-1, 1))
    
    p3 <- ggplot(res_bias_lfe[[3]][[2]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Years of life with dementia (PBS 3-6 years)",
           x = "Model",
           y = "Rel Bias") +
      theme_bw() +
      geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 9)) +
      scale_y_continuous(limits = c(-1, 1))
    
    p4 <- ggplot(res_bias_lfe[[4]][[2]], aes(x = model, y = lfe, color = model)) +
      geom_point(size = 4) +
      labs(title = "Years of life with dementia (EHR)",
           x = "Model",
           y = "Rel Bias") +
      theme_bw() +
      geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 9)) +
      scale_y_continuous(limits = c(-1, 1))
    
    (p1 + p2 + p3 + p4) +
      plot_layout(ncol = 2, nrow = 2, guides = "collect") +
      theme(plot.margin = margin(10, 10, 10, 10),
            legend.position = "right")
  }
}
