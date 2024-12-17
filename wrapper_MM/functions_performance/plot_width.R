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

  ggplot(df_long, aes(x = model, y = se_value, color = model)) +
    geom_point(size=3, shape= "diamond") +
    facet_grid(covariate~transition, scales = "fixed", labeller = as_labeller(c(
      "1" = "Dementia-free -> Dementia",  
      "2" = "Dementia-free -> Death",
      "3" = "Dementia -> Death",
      "cov1"="beta1",
      "cov2"="beta2",
      "cov3"="beta3"
    ))) +    
    #coord_cartesian(ylim = c(0, max(se_value)))+
    labs(title = titles[scheme - 1],  
         x = "Model",
         y = "Mean SE") +
    theme_bw() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme(axis.text.x = element_blank(), 
          legend.position = "right",
          axis.title.x = element_blank()) +
    scale_color_viridis_d() 
  
 
  
}