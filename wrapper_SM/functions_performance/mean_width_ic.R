mean_width_ic <- function(data, scheme){
  temp <- data[[scheme-1]]
  temp <- as.data.frame(temp)
  desired_order <- c("flexsurv_EO", "coxph", "flexsurv","imputation")
  temp$model <- factor(temp$model, levels = desired_order, ordered = TRUE)
  temp <- temp[order(temp$model), ]
  temp <- temp %>%
    mutate(across(1:3, as.numeric))
  
  mean_width <- temp %>%
    group_by(model, covariate, ic) %>%
    summarise(
      across(c(1, 2, 3), 
             ~ round(mean(.x, na.rm = TRUE), 3)), 
      .groups = 'drop'
    )
  
  return (mean_width)
}