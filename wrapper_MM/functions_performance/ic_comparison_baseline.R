ic_comparison_baseline <- function(data, scheme) {
  temp <- data[[scheme-1]]
  temp <- as.data.frame(temp)
  
  desired_order <- c("flexsurv_EO","flexsurv", "msm_age", "nhm", "imputation")
  temp$model <- factor(temp$model, levels = desired_order, ordered = TRUE)
  temp <- temp[order(temp$model), ]
  
  
  temp <- temp %>%
    mutate(across(1:2, as.numeric))
  
  #ci using
  ci_baseline<- temp %>%
    group_by(model, transition) %>%
    summarise(
      across(c(rate, shape), 
             list(
               mean = ~ round(mean(.x, na.rm = TRUE), 5),
               lower_ci = ~ round(mean(.x, na.rm = TRUE) - qt(0.975, df = sum(!is.na(.x)) - 1) * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))), 5),
               upper_ci = ~ round(mean(.x, na.rm = TRUE) + qt(0.975, df = sum(!is.na(.x)) - 1) * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))), 5)
             )
      ),
      .groups = 'drop'
    )
  
  return(ci_baseline)
}
