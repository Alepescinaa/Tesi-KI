compute_CI <- function(T1,W,B,m,theta_bar){
  # t-student
  # SE <- sqrt(diag(T1))
  # nu <- numeric()
  # for (i in 1:5){
  #   nu[i] <- (m - 1) * (1 + W[i,i] / B[i,i])^2 / (1 + (1 / (m + 1)) * (W[i,i] / B[i,i])^2)
  # }
  # alpha <- 0.05
  # t_value <- numeric()
  # for (i in 1:5){
  #   t_value[i] <- qt((1 - alpha / 2), df = nu[i])  
  # }
  # 
  # lower_bound <- numeric()
  # upper_bound <- numeric()
  # conf_int <- matrix(0,5,2)
  # for (i in 1:5){
  #   lower_bound[i]<- theta_bar[i] - t_value[i] * SE[i]
  #   upper_bound[i] <- theta_bar[i] + t_value[i] * SE[i]
  #   conf_int[i,] <- cbind(lower_bound[i], upper_bound[i])
  # }
  # colnames(conf_int) <- c("L95%","U95%")
  
  # normal 
  SE <- sqrt(diag(T1)) 
  alpha <- 0.05
  
  z_value <- qnorm(1 - alpha / 2) 
  
  lower_bound <- theta_bar - z_value * SE
  upper_bound <- theta_bar + z_value * SE
  
  conf_int <- cbind(lower_bound, upper_bound)
  colnames(conf_int) <- c("L95%", "U95%")
  return(conf_int)
}