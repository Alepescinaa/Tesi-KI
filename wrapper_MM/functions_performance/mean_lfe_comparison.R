mean_lfe_comparison <- function(data, scheme){
  temp <- data[[scheme-1]]
  temp <- as.data.frame(temp)
  desired_order <- c("flexsurv_EO", "coxph", "flexsurv", "msm", "msm_age", "nhm", "imputation")
  temp$model <- factor(temp$model, levels = desired_order, ordered = TRUE)
  temp <- temp[order(temp$model), ]
  temp$lfe <- as.numeric(temp$lfe)
  d1 <- temp %>%
    group_by(model,seed) %>%
    slice(1) %>%
    ungroup() 
  d2 <- temp %>%
    group_by(model,seed) %>%
    slice(2) %>%
    ungroup() 
  
  tls_Health <- d1 %>%
    group_by(model) %>%
    summarise(lfe = mean(lfe, na.rm = TRUE)) 
  
  tls_Dem <- d2 %>%
    group_by(model) %>%
    summarise(lfe = mean(lfe, na.rm = TRUE)) 
  
  return (list(tls_Health,tls_Dem))
}