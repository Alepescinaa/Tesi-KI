mean_bias_comparison <- function(bias_all_schemes, scheme){
  temp <- bias_all_schemes[[scheme-1]]
  temp <- as.data.frame(temp)
  desired_order <- c("flexsurv_EO", "coxph", "flexsurv", "imputation")
  temp$model <- factor(temp$model, levels = desired_order, ordered = TRUE)
  temp <- temp[order(temp$model), ]
  temp <- temp %>%
    mutate(across(1:8, as.numeric))
  
  mean_bias <- temp %>%
    group_by(model, transition) %>%
    summarise(
      across(c(rate, shape, cov1, cov2, cov3, `exp(cov1)`, `exp(cov2)`, `exp(cov3)`), 
             ~ round(mean(.x, na.rm = TRUE), 3)), 
      .groups = 'drop'
    )
  
  return (mean_bias)
}