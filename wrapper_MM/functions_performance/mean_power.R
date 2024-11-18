mean_power <- function(results, scheme){
  temp <- results
  temp <- as.data.frame(temp)
  desired_order <- c("flexsurv_EO", "coxph", "flexsurv", "msm", "msm_age", "nhm", "imputation")
  temp$model <- factor(temp$model, levels = desired_order, ordered = TRUE)
  temp <- temp[order(temp$model), ]
  temp <- temp %>%
    mutate(across(c(1:3,6), as.numeric))
  
  mean_bias <- temp %>%
    group_by(model, transition) %>%
    summarise(
      across(c(cov1, cov2, cov3),
             ~ round(mean(.x, na.rm = TRUE)*100, 3)), 
      .groups = 'drop'
    )
  
  return (mean_bias)
}