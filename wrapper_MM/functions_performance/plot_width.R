plot_se <- function (data,scheme, titles){
  df <- data[[scheme-1]]
  df <- as.data.frame(df)
  
  df <- df %>%
    # mutate( beta1=cov1, beta2=cov2, beta3=cov3, cov1=NULL, cov2=NULL, cov3=NULL) %>%
    mutate(
      model = case_when(
        model == "flexsurv_EO" ~ "0",
        model == "coxph" ~ "a",
        model == "flexsurv" ~ "b",
        model == "msm" ~ "c",
        model == "msm_age" ~ "d",
        model == "nhm" ~ "e",
        model == "imputation" ~ "f",
        TRUE ~ model # Keep other values unchanged if they exist
      )
    )
  
  
  df_long <- df %>%
    pivot_longer(cols = c("1", "2", "3"),
                 names_to = "transition",
                 values_to = "se_value")
  
  latex_labels <- c(
    "1" = "Transition 1",
    "2" = "Transition 2",
    "3" = "Transition 3",
    "cov1" = "beta[1]",
    "cov2" = "beta[2]",
    "cov3" = "beta[3]"
  )

  ggplot(df_long, aes(x = model, y = se_value, color = model)) +
    geom_point(size=3, shape= "diamond") +
    facet_grid(covariate~transition, scales = "fixed", labeller = labeller(
      covariate = as_labeller(latex_labels, label_parsed),
      transition = as_labeller(latex_labels)
    )) +    
    #coord_cartesian(ylim = c(0, max(se_value)))+
    labs(title = titles[scheme - 1],  
         x = "Model",
         y = "Standard Error") +
    theme_bw() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_viridis_d() 
  
 
  
}