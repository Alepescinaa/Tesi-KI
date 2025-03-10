mean_standard_error <- function(data, scheme){
  temp <- data[[scheme-1]]
  temp <- as.data.frame(temp)
  desired_order <- c("flexsurv_EO", "coxph", "flexsurv", "msm", "msm_age", "nhm", "imputation")
  temp$model <- factor(temp$model, levels = desired_order, ordered = TRUE)
  temp <- temp[order(temp$model), ]
  temp <- temp %>%
    mutate(across(1:3, as.numeric))
  
  se <- temp %>%
    group_by(model, covariate) %>%
    summarise(
      across(c(1, 2, 3), 
             ~ round(mean(.x, na.rm = TRUE), 3)), 
      .groups = 'drop'
    )
  
  return (se)
}