plot_bias <- function (data, scheme, titles){
  
  df <- data[[scheme-1]]
  
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
  
  
  df_long_mean<- df %>%
    pivot_longer(cols = c("rate_mean", "shape_mean", "cov1_mean", "cov2_mean", "cov3_mean"),
                 names_to = "parameter",
                 values_to = "bias_mean") %>% 
    select(model,transition,parameter,bias_mean) %>% 
    mutate(parameter=gsub("_mean","",parameter))
  
  df_long_upper<- df %>%
    pivot_longer(cols = c("rate_upper_ci", "shape_upper_ci", "cov1_upper_ci", "cov2_upper_ci", "cov3_upper_ci"),
                 names_to = "parameter",
                 values_to = "bias_upper")%>% 
    select(model,transition,parameter,bias_upper)%>% 
    mutate(parameter=gsub("_upper_ci","",parameter))
  
  df_long_lower<- df %>%
    pivot_longer(cols = c("rate_lower_ci", "shape_lower_ci", "cov1_lower_ci", "cov2_lower_ci", "cov3_lower_ci"),
                 names_to = "parameter",
                 values_to = "bias_lower")%>% 
    select(model,transition,parameter,bias_lower)%>% 
    mutate(parameter=gsub("_lower_ci","",parameter))
  
  df_long<-df_long_mean %>% left_join(df_long_lower) %>% left_join(df_long_upper)
  
  latex_labels <- c(
    "1" = "Transition 1",
    "2" = "Transition 2",
    "3" = "Transition 3",
    "cov1" = "beta[1]",
    "cov2" = "beta[2]",
    "cov3" = "beta[3]"
  )

  cov_plot <- ggplot(data = df_long %>% 
                filter(parameter %in% c( "cov1", "cov2", "cov3")),
              aes(x = model, y = bias_mean, color = model)) +
    geom_point(size = 3) +
    facet_grid(parameter~transition, scales = "free",  labeller = labeller(
      parameter = as_labeller(latex_labels, label_parsed),
      transition = as_labeller(latex_labels)
    )) +    
    labs(title = titles[scheme - 1],  
         x = "Model",
         y = "Bias") +
    geom_errorbar(aes(x=model, ymin= bias_lower, ymax= bias_upper, color= model, width= 0.3))+
    theme_bw() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_line(data = df_long %>% 
                filter(parameter %in% c( "cov1", "cov2", "cov3")),
                aes(group = transition), 
                size = 0.5, 
                linetype = "dashed", 
                color = "grey", 
                na.rm = TRUE)+ 
    scale_color_viridis_d() 
  
  latex_labels <- c(
    "1" = "Transition 1",
    "2" = "Transition 2",
    "3" = "Transition 3",
    "shape" = "mu",
    "rate" = "lambda"
  )
  
  baseline_plot <- ggplot(data = df_long %>% 
                filter(parameter %in% c( "rate", "shape"),
                       model %in% c( "0", "b", "d", "e", "f")),
              aes(x = model, y = bias_mean, color = model)) +
    geom_point(size = 3) +
    facet_grid(parameter~transition, scales = "free", labeller = labeller(
      parameter = as_labeller(latex_labels, label_parsed),
      transition = as_labeller(latex_labels)
    )) +    
    labs(title = titles[scheme - 1],  
         x = "Model",
         y = "Bias") +
    theme_bw() +
    geom_errorbar(aes(x=model, ymin=bias_lower, ymax= bias_upper, color= model, width=0.3))+
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_line(data = df_long %>% 
                filter(parameter %in% c( "rate", "shape"),
                       model %in% c( "0", "b", "d", "e", "f")),
              aes(group = transition), 
              size = 0.5, 
              linetype = "dashed", 
              color = "grey", 
              na.rm = TRUE) +
    scale_color_manual(values = c(
      "0" = viridis::viridis(7)[1],  
      "b" = viridis::viridis(7)[3],  
      "d" = viridis::viridis(7)[5],  
      "e" = viridis::viridis(7)[6],  
      "f" = viridis::viridis(7)[7]    
    ))
  
  return(list(cov_plot,baseline_plot))

  
}
