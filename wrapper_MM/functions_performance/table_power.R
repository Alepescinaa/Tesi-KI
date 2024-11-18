table_power <- function(significant_covs, power, scheme){
  df_merged <- power[[scheme-1]] %>%
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
    column_spec(3, background = ifelse(df_merged$cov1_color != "", "yellow", "white")) %>%
    column_spec(4, background = ifelse(df_merged$cov2_color != "", "yellow", "white")) %>%
    column_spec(5, background = ifelse(df_merged$cov3_color != "", "yellow", "white"))
  
}