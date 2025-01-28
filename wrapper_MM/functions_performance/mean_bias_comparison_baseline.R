mean_bias_comparison_baseline <- function(data, scheme){
  temp <- data[[scheme-1]]
  temp <- as.data.frame(temp)
  desired_order <- c("flexsurv_EO", "flexsurv",  "msm_age", "nhm", "imputation")
  temp$model <- factor(temp$model, levels = desired_order, ordered = TRUE)
  temp <- temp[order(temp$model), ]
  temp <- temp %>%
    mutate(across(1:2, as.numeric))
  
  mean_bias <- temp %>%
    group_by(model, transition) %>%
    summarise(
      across(c(rate, shape), 
             ~ round(mean(.x, na.rm = TRUE), 3)), 
      .groups = 'drop'
    )
  
  return (mean_bias)
}