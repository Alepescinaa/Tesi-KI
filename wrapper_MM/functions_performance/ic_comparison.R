ic_comparison <- function(data, scheme) {
  temp <- data[[scheme-1]]
  temp <- as.data.frame(temp)
  
  desired_order <- c("flexsurv_EO", "coxph", "flexsurv", "msm", "msm_age", "nhm", "imputation")
  temp$model <- factor(temp$model, levels = desired_order, ordered = TRUE)
  temp <- temp[order(temp$model), ]
  

  temp <- temp %>%
    mutate(across(1:5, as.numeric))
  
  #ci using
  ci_cov <- temp %>%
    group_by(model, transition) %>%
    summarise(
      across(c(rate, shape, cov1, cov2, cov3), 
             list(
               mean = ~ round(mean(.x, na.rm = TRUE), 5),
               lower_ci = ~ round(mean(.x, na.rm = TRUE) - qt(0.975, df = sum(!is.na(.x)) - 1) * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))), 5),
               upper_ci = ~ round(mean(.x, na.rm = TRUE) + qt(0.975, df = sum(!is.na(.x)) - 1) * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))), 5)
             )
      ),
      .groups = 'drop'
    )
  
  return(ci_cov)
}
