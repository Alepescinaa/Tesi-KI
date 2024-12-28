
plot_lfe_bias <- function(ci_lfe) {
    
    df_all <- rbind(
      transform(ci_lfe[[1]], scenario = "1 year observation scheme"),
      transform(ci_lfe[[2]], scenario = "3 years observation scheme"),
      transform(ci_lfe[[3]], scenario = "3-6 years observation scheme"),
      transform(ci_lfe[[4]], scenario = "irregular observation scheme")
    )
    
    df_all <- df_all %>%
      mutate(
        model = case_when(
          model == "flexsurv_EO" ~ "0",
          model == "coxph" ~ "a",
          model == "flexsurv" ~ "b",
          model == "msm" ~ "c",
          model == "msm_age" ~ "d",
          model == "nhm" ~ "e",
          model == "imputation" ~ "f",
          TRUE ~ model 
        )
      )
    

    ggplot(df_all, aes(x = model, y = lfe_mean, fill = as.factor(dem))) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.4) +
      geom_errorbar(aes(ymin = lfe_lower_ci, ymax = lfe_upper_ci),
                    position = position_dodge(width = 0.7), width = 0.25) +
      labs(
        title = "Relative Bias of Time Spent in State",
        x = NULL,
        y = "Relative Bias (%)",
        fill = NULL
      ) +
      scale_y_continuous(limits = c(-1, 0.7)) +
      scale_fill_manual(values = c("0" = "skyblue", "1" = "purple"),
                        labels = c("Dementia-free", "Dementia")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.position = "top", 
        legend.margin = margin(t = 10),  
        panel.spacing = unit(1, "lines") 
      ) +
      facet_wrap(~scenario, ncol = 2)
    
}
    
