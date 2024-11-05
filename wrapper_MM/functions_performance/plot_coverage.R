plot_coverage <- function (scheme, titles){
  df <- res_cov[[scheme-1]]
  df <- as.data.frame(df)
  desired_order <- c("flexsurv_EO", "coxph", "flexsurv", "msm", "msm_age", "nhm", "imputation")
  df$model <- factor(df$model, levels = desired_order, ordered = TRUE)
  df <- df[order(df$model), ]
  df_long <- df %>%
    pivot_longer(cols = c("rate", "shape", "cov1", "cov2", "cov3"),
                 names_to = "parameter",
                 values_to = "coverage")
  
  ggplot(df_long, aes(x = model, y = coverage, color = model, shape = transition)) +
    geom_point(size = 3) +
    facet_wrap(~ parameter, scales = "fixed") +  
    labs(title = titles[scheme - 1],  
         x = "Model",
         y = "Coverage") +
    theme_bw() +
    geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", size = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")  
  

}