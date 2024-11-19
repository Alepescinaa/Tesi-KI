mean_coverage_comparison <- function(coverage_all_schemes, scheme){
  
  temp <- coverage_all_schemes[[scheme-1]]
  temp <- as.data.frame(temp)
  temp <- temp %>%
    mutate(across(1:5, as.numeric))
  
  mean_coverage <- temp %>%
    group_by(model, transition) %>%
    summarise(
      across(c(rate, shape, cov1, cov2, cov3), 
             ~ round(mean(.x, na.rm = TRUE), 3)), 
      .groups = 'drop'
    )
  
  return(mean_coverage)
}


