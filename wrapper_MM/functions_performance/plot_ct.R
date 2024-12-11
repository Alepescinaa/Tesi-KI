plot_ct <- function(scheme, titles, convergence){
  temp <- ct_all_schemes[[scheme-1]]
  model_names <- c("coxph_ct", "flexsurv_ct", "msm_ct", "msm_age_ct", "nhm_ct","imputation_ct")
  colnames(temp) <- model_names
  temp <- cbind(temp,convergence)
  temp$coxph_ct <- ifelse(temp$coxph == 2, temp$coxph_ct, NA)
  temp$flexsurv_ct <- ifelse(temp$flexsurv == 2, temp$flexsurv_ct, NA)
  temp$msm_ct <- ifelse(temp$msm == 2, temp$msm_ct, NA)
  temp$msm__age_ct <- ifelse(temp$msm_age == 2, temp$msm_age_ct, NA)
  temp$nhm_ct <- ifelse(temp$nhm == 2, temp$nhm_ct, NA)
  temp$imputation_ct <- ifelse(temp$imputation == 2, temp$imputation_ct, NA)
  mean_times <- colMeans(temp[1:6], na.rm=T)
  mean_times_df <- data.frame(model = names(mean_times), time = mean_times)
  ggplot(mean_times_df, aes(x = model, y = time, fill = model)) +
    geom_bar(stat = "identity", color = "black") +
    labs(title = titles[scheme-1],
         x = "Model",
         y = "Mean Computational Time (sec)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(limits = c(0, 1500))
  
  
}
