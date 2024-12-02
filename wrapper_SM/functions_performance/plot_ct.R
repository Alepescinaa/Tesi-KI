plot_ct <- function(scheme, titles, convergence){
  temp <- ct_all_schemes[[scheme-1]]
  temp <- as.data.frame(temp)
  temp <- temp %>%
    mutate(V3 = NULL)
  model_names <- c("coxph_ct", "flexsurv_ct","imputation_ct")
  colnames(temp) <- model_names
  temp <- cbind(temp,convergence)
  temp$coxph_ct <- ifelse(temp$coxph == 2, temp$coxph_ct, NA)
  temp$flexsurv_ct <- ifelse(temp$flexsurv == 2, temp$flexsurv_ct, NA)
  temp$imputation_ct <- ifelse(temp$imputation == 2, temp$imputation_ct, NA)
  mean_times <- colMeans(temp[1:3], na.rm=T)
  mean_times_df <- data.frame(model = names(mean_times), time = mean_times)
  ggplot(mean_times_df, aes(x = model, y = time, fill = model)) +
    geom_bar(stat = "identity", color = "black") +
    labs(title = titles[scheme-1],
         x = "Model",
         y = "Mean Computational Time (sec)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(limits = c(0, 30))
  print(mean_times)
  
  
}
