plot_distribution <- function (data, scheme){
  
  df <- data[[scheme-1]]
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
    ) %>%
    select(model, transition, rate, shape, cov1, cov2, cov3) %>%
    mutate( beta1=cov1, beta2=cov2, beta3=cov3, cov1=NULL, cov2=NULL, cov3=NULL)
  
  df_long<- df %>%
    pivot_longer(cols = c("rate", "shape", "beta1", "beta2", "beta3"),
                 names_to = "parameter",
                 values_to = "estimates") 
  df_long$estimates <- as.numeric( df_long$estimates)
  
  df_long1 <- df_long[df_long$transition==1,]
  df_long2 <- df_long[df_long$transition==2,]
  df_long3 <- df_long[df_long$transition==3,]

  p1 <- ggplot(df_long1 %>% 
           filter(parameter %in% c( "beta1", "beta2", "beta3")), aes(x = estimates, y = model, fill = model)) +
    geom_density_ridges() +
    theme_ridges() + 
    facet_grid(parameter  ~ ., scales = "free")+
    # geom_violin(trim = TRUE) +  
    labs(
      x =NULL,
      y=NULL,
      title = "Distribution of estimates for transition Dementia-free -> Dementia"
    ) +
    scale_fill_manual(values = c(
      "0" = viridis::viridis(7)[1], 
      "a" = viridis::viridis(7)[2],
      "b" = viridis::viridis(7)[3], 
      "c" = viridis::viridis(7)[4],
      "d" = viridis::viridis(7)[5],  
      "e" = viridis::viridis(7)[6],  
      "f" = viridis::viridis(7)[7]    
    ))+
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  p2 <- ggplot(df_long2 %>% 
           filter(parameter %in% c( "beta1", "beta2", "beta3")), aes(x = estimates, y = model, fill = model)) +
    geom_density_ridges() +
    theme_ridges() + 
    facet_grid(parameter  ~ ., scales = "free")+
    # geom_violin(trim = TRUE) +  
    labs(
      x =NULL,
      y=NULL,
      title = "Distribution of estimates for transition Dementia-free -> Death"
    ) +
    scale_fill_manual(values = c(
      "0" = viridis::viridis(7)[1], 
      "a" = viridis::viridis(7)[2],
      "b" = viridis::viridis(7)[3], 
      "c" = viridis::viridis(7)[4],
      "d" = viridis::viridis(7)[5],  
      "e" = viridis::viridis(7)[6],  
      "f" = viridis::viridis(7)[7]    
    ))+
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  p3 <- ggplot(df_long3 %>% 
           filter(parameter %in% c( "beta1", "beta2", "beta3")), aes(x = estimates, y = model, fill = model)) +
    geom_density_ridges() +
    theme_ridges() + 
    facet_grid(parameter  ~ ., scales = "free")+
    # geom_violin(trim = TRUE) +  
    labs(
      x =NULL,
      y=NULL,
      title = "Distribution of estimates for transition Dementia-> Death"
    ) +
    scale_fill_manual(values = c(
      "0" = viridis::viridis(7)[1], 
      "a" = viridis::viridis(7)[2],
      "b" = viridis::viridis(7)[3], 
      "c" = viridis::viridis(7)[4],
      "d" = viridis::viridis(7)[5],  
      "e" = viridis::viridis(7)[6],  
      "f" = viridis::viridis(7)[7]    
    ))+
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
 
  return(list(p1,p2,p3))
  
  
}
