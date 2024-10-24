mean_bias_comparison <- function(bias_all_schemes, scheme){
  temp <- bias_all_schemes[[scheme-1]]
  temp <- as.data.frame(temp)
  temp <- temp %>%
    mutate(across(1:8, as.numeric))
  
  nhm_missing <- temp %>%
    filter(model == "nhm" & 
             is.na(rate) & 
             is.na(shape) & 
             is.na(cov1) & 
             is.na(cov2) & 
             is.na(cov3) & 
             is.na(`exp(cov1)`) & 
             is.na(`exp(cov2)`) & 
             is.na(`exp(cov3)`)
    ) %>%
    select(seed) %>%
    distinct()
  nhm_missing <- as.numeric(nhm_missing$seed)
  
  mean_bias <- temp %>%
    group_by(model, transition) %>%
    summarise(
      across(c(rate, shape, cov1, cov2, cov3, `exp(cov1)`, `exp(cov2)`, `exp(cov3)`), 
             ~ round(mean(.x, na.rm = TRUE), 3)), 
      .groups = 'drop'
    )
  
  return (list(mean_bias=mean_bias, nhm_missing=nhm_missing))
}