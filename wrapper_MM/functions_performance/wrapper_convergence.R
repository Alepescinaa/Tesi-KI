wrapper_convergence <- function(n_pats, scheme, seed){
  temp <- vector(mode = "list", length = 100)
  res_checking <- vector(mode = "list", length = 100)
  res_hessian <- vector(mode = "list", length = 100)

  for (seed in 1:100){
    setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI/wrapper_MM")
    temp[[seed]] <- check_convergence(n_pats, scheme, seed)
    res_checking[[seed]] <- temp[[seed]][[1]]
    res_hessian[[seed]] <- temp[[seed]][[2]]
  }
  
  res_checking <-  do.call(rbind, res_checking)
  res_hessian <-do.call(rbind, res_hessian)
  res_checking <- as.data.frame(res_checking)
  res_hessian <- as.data.frame(res_hessian)
  res_checking <- res_checking %>%
    mutate(across(1:6, as.numeric))
  res_checking <- rbind(res_checking, colMeans(res_checking)) 
  res_hessian <- res_hessian %>%
    mutate(across(1:6, as.numeric))
  res_hessian <- rbind(res_hessian, colMeans(res_hessian)) 
  return(list(res_checking = res_checking, res_hessian = res_hessian))
}