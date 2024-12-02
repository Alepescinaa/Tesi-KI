plot_coverage <- function (scheme, titles){
  df <- res_cov_ic[[scheme-1]]
  df <- as.data.frame(df)
  desired_order <- c("flexsurv_EO", "coxph", "flexsurv", "imputation")
  df$model <- factor(df$model, levels = desired_order, ordered = TRUE)
  df <- df[order(df$model), ]
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
  
  ggplot(df_long, aes(x = model, y = coverage_mean, color = model, shape = transition)) +
    geom_point(size = 3) +
    geom_errorbar(aes(x=model, ymin= coverage_lower, ymax= coverage_upper, color= model))+
    facet_grid(transition~ parameter, scales = "free") +  
    labs(title = titles[scheme - 1],  
         x = "Model",
         y = "Coverage") +
    theme_bw() +
    geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", size = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")  
  

}