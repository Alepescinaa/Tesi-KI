prepare_data_boxplot <- function(scheme){
  temp <- bias_all_schemes[[scheme-1]]
  temp <- as.data.frame(temp)
  desired_order <- c("flexsurv_EO", "coxph", "flexsurv", "msm", "msm_age", "nhm", "imputation")
  temp$model <- factor(temp$model, levels = desired_order, ordered = TRUE)
  temp <- temp[order(temp$model), ]
  temp1 <- temp[temp$transition==1,]
  temp2 <- temp[temp$transition==2,]
  temp3 <- temp[temp$transition==3,]
  
  df_long1 <- temp1 %>%
    pivot_longer(cols = c("rate", "shape", "cov1", "cov2", "cov3"),
                 names_to = "parameter",
                 values_to = "bias")
  df_long1 <- df_long1 %>%
    select(-c('exp(cov1)', 'exp(cov2)', 'exp(cov3)'))
  df_long1$bias <- as.numeric(df_long1$bias)
  
  df_long2 <- temp2 %>%
    pivot_longer(cols = c("rate", "shape", "cov1", "cov2", "cov3"),
                 names_to = "parameter",
                 values_to = "bias")
  df_long2 <- df_long2 %>%
    select(-c('exp(cov1)', 'exp(cov2)', 'exp(cov3)'))
  df_long2$bias <- as.numeric(df_long2$bias)
  

  df_long3 <- temp3 %>%
    pivot_longer(cols = c("rate", "shape", "cov1", "cov2", "cov3"),
                 names_to = "parameter",
                 values_to = "bias")
  df_long3 <- df_long3 %>%
    select(-c('exp(cov1)', 'exp(cov2)', 'exp(cov3)'))
  df_long3$bias <- as.numeric(df_long3$bias)
  
  return(list(df_long1,df_long2,df_long3))
}