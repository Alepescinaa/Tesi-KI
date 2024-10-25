wrapper_convergence <- function(n_pats, scheme, seed){
  res_checking <- vector(mode = "list", length = 100)
  for (seed in 1:100){
    setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI/wrapper_MM")
    res_checking[[seed]] <- check_convergence(n_pats, scheme, seed)
  }
  
  res_checking <-  do.call(rbind, res_checking)
  res_checking <- as.data.frame(res_checking)
  res_checking <- res_checking %>%
    mutate(across(1:4, as.numeric))
  res_checking <- rbind(res_checking, colMeans(res_checking)) 
  return(res_checking)
}