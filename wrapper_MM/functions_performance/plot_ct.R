plot_ct <- function(scheme, titles){
  temp <- ct_all_schemes[[scheme-1]]
  model_names <- c("coxph", "flexsurv", "msm", "msm_age", "nhm","imputation")
  colnames(temp) <- model_names
  temp <- cbind(temp,convergence_schemes[[scheme-1]][1:100,])
  temp$msm <- ifelse(temp$converged_msm == 0, NA, temp$msm)
  temp$msm_age<- ifelse(temp$converged_msm_age == 0, NA, temp$msm_age)
  temp$nhm<- ifelse(temp$nhm_computed == 0, NA, temp$nhm)
  temp$coxph<- ifelse(temp$converged_coxph == 0, NA, temp$coxph)
  mean_times <- colMeans(temp[1:6], na.rm=T)
  mean_times_df <- data.frame(model = names(mean_times), time = mean_times)
  ggplot(mean_times_df, aes(x = model, y = time, fill = model)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = c("red", "blue", "green", "purple", "black", "orange")) +
    labs(title = titles[scheme-1],
         x = "Model",
         y = "Mean Computational Time (sec)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  
}