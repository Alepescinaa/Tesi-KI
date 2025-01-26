plot_coverage <- function (scheme, titles){
  df <- res_cov_ic[[scheme-1]]
  df <- df %>%
    mutate(across(c(3:17), ~ pmax(pmin(., 1), 0)))
  df <- as.data.frame(df)
  df <- df %>%
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
  
  df <- df %>%
    mutate(across(
      .cols = !c(model, transition),
      .fns = ~ . * 100             
    ))
  
  df_long_mean<- df %>%
    pivot_longer(cols = c("rate_mean", "shape_mean", "cov1_mean", "cov2_mean", "cov3_mean"),
                 names_to = "parameter",
                 values_to = "coverage_mean") %>% 
    select(model,transition,parameter,coverage_mean) %>% 
    mutate(parameter=gsub("_mean","",parameter))
  
  df_long_upper<- df %>%
    pivot_longer(cols = c("rate_upper_ci", "shape_upper_ci", "cov1_upper_ci", "cov2_upper_ci", "cov3_upper_ci"),
                 names_to = "parameter",
                 values_to = "coverage_upper")%>% 
    select(model,transition,parameter,coverage_upper)%>% 
    mutate(parameter=gsub("_upper_ci","",parameter))
  
  df_long_lower<- df %>%
    pivot_longer(cols = c("rate_lower_ci", "shape_lower_ci", "cov1_lower_ci", "cov2_lower_ci", "cov3_lower_ci"),
                 names_to = "parameter",
                 values_to = "coverage_lower")%>% 
    select(model,transition,parameter,coverage_lower)%>% 
    mutate(parameter=gsub("_lower_ci","",parameter))
  
  df_long<-df_long_mean %>% left_join(df_long_lower) %>% left_join(df_long_upper)
  
  latex_labels <- c(
    "1" = "Transition 1",
    "2" = "Transition 2",
    "3" = "Transition 3",
    "cov1" = "beta[1]",
    "cov2" = "beta[2]",
    "cov3" = "beta[3]",
    "rate" ="lambda",
    "shape"= "gamma" 
  )
  
  
  ggplot(df_long, aes(x = model, y = coverage_mean, color = model)) +
    geom_point(size = 3) +
    geom_errorbar(aes(x=model, ymin= coverage_lower, ymax= coverage_upper, color= model))+
    facet_grid(transition~ parameter, scales = "fixed",   labeller = labeller(
      parameter = as_labeller(latex_labels, label_parsed),
      transition = as_labeller(latex_labels)
    )) +    
    labs(title = titles[scheme - 1],  
         x = "Model",
         y = "Coverage (%)") +
    theme_bw() +
    geom_hline(yintercept = 95, color = "red", linetype = "dashed", size = 0.5) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_viridis_d()  
  

}