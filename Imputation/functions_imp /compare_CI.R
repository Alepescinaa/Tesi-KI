compare_CI <- function(data1,data2){
  df_imp <- as.data.frame(data1)
  df_EO <- as.data.frame(data2)
  
  df_imp$Parameter <- paste( c("shape", "rate" , "cov1" , "cov2","cov3") )
  df_EO$Parameter <-paste( c("shape", "rate" , "cov1" , "cov2","cov3") )
  
  combined_df <- df_imp %>%
    rename(LCI = "L95%", UCI = "U95%") %>%
    full_join(df_EO %>% rename(LCI = "L95%", UCI = "U95%"), by = "Parameter", suffix = c("_A", "_B"))
  
  create_plot <- function(row){
    p <- ggplot(row, aes(x = Parameter)) +
      geom_segment(aes(xend = Parameter, y = LCI_A, yend = UCI_A), color = "blue", size = 1) +
      geom_segment(aes(xend = Parameter, y = LCI_B, yend = UCI_B), color = "red", size = 1, linetype = "dashed") + 
      labs(title = "Confidence Intervals for Parameters",
           x = "Parameters",
           y = "Confidence Interval") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  lapply(1:nrow(combined_df), function(i) create_plot(combined_df[i, ]))
  
}