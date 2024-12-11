table_power <- function(significant_covs, data, scheme){
  data <- data[[scheme-1]]
  # data <- data %>%
  #   mutate(
  #     cov1 = (100 - cov1),
  #     cov2 = (100 - cov2),
  #     cov3 = (100 - cov3)
  #   )
  # 
  max_power <-  data %>%
    pivot_longer(cols = starts_with("cov"), 
                 names_to = "covariate", 
                 values_to = "value") %>% 
    group_by(transition, covariate) %>%
    slice_max(value, with_ties = FALSE) %>% 
    ungroup()
  max_power$value[max_power$transition==1 & max_power$covariate=="cov1"] <- NA
  max_power$value[max_power$transition==3 & max_power$covariate=="cov2"] <- NA
  
  print(max_power)

  
  df_merged <- data %>%
    left_join(significant_covs, by = c( "transition"))
  
  df_merged <- df_merged %>%
    mutate(
      cov1_color = ifelse(cov1.y == 1, "background-color: yellow;", ""),
      cov2_color = ifelse(cov2.y == 1, "background-color: yellow;", ""),
      cov3_color = ifelse(cov3.y == 1, "background-color: yellow;", "")
    )
  
  df_merged %>%
    select(model, transition, cov1.x, cov2.x, cov3.x) %>%
    kable("html", escape = FALSE) %>%
    kable_styling() %>%
    column_spec(3, background = ifelse(df_merged$cov1_color != "", "yellow", "lightgreen")) %>%
    column_spec(4, background = ifelse(df_merged$cov2_color != "", "yellow", "lightgreen")) %>%
    column_spec(5, background = ifelse(df_merged$cov3_color != "", "yellow", "lightgreen")) %>%
    add_header_above(c(" "=2, "Power vs Error type I" = 3))
  

  
}