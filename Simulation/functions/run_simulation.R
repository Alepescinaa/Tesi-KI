run_simulation <- function(n_pats, fits_wei, model_type, seed, followup, meanlog, sdlog ) {
  set.seed(seed)
  res <- simulation(n_pats, fits_wei, model_type, meanlog, sdlog) 
  sim_data <- res[[1]]
  extremes<- res[[2]]
  temp <- sim_data$sim_disease(max_t=120, max_age=120) 
  temp[, c("sample", "strategy_id", "grp_id") := NULL]
  
  delete <- which(temp$to==2)
  temp <- temp[-delete,]
  temp$from <- temp$from-1
  temp$to <- temp$to-1
  
  cov_to_add <- as.data.frame(sim_data$input_data$X[[2]][[2]][,-1])
  #cov_to_add[,1] <- cov_to_add[,1]*(extremes[2]-extremes[1])+extremes[1]
  cov_to_add$patient_id <- 1:n_pats
  temp <- temp %>%
    left_join(cov_to_add, by = "patient_id")
  
  
  
  # t <- seq(60,100,2)
  # probs <- sim_data$sim_stateprobs(t)
  # probs[, c("sample", "strategy_id", "grp_id") := NULL]
  
  
  keep_index <- temp$patient_id[temp$from==1 & temp$to==2]
  
  temp <- apply_right_censoring(temp)
  list_temp <- replicate(5, temp, simplify = FALSE)
  list_temp[[1]]$age <- list_temp[[1]]$time_start
  list_temp[[1]] <- list_temp[[1]] %>%
    group_by(patient_id) %>%
    mutate(age = ifelse(n() == 2, first(time_start), age)) %>%
    ungroup()
  list_temp[[1]]$onset <- ifelse((list_temp[[1]]$from==2 & list_temp[[1]]$to==3),1,0) 
  list_temp[[1]]$onset_age <- list_temp[[1]]$death_time
  list_temp[[1]]$onset_age[list_temp[[1]]$from==2 & list_temp[[1]]$to==3] <- list_temp[[1]]$time_stop[list_temp[[1]]$from==1 & list_temp[[1]]$to==2]
  list_temp[[1]] <- list_temp[[1]] %>%
    group_by(patient_id) %>%
    slice(ifelse(n() > 1, 2, 1)) %>%
    ungroup()
  list_temp[[1]] <- list_temp[[1]][,c(1,7:9,10:15)]
  
  apply_scheme <- function(j){
    res <- list()
    for (i in 1:n_pats ){
      temp_id <- list_temp[[j+1]][list_temp[[j+1]]$patient_id == i, ]
      res[[i]] <- interval_censoring_id(temp_id, i, scheme = j, followup, keep_index)
    }
    list_temp[[j+1]] <- do.call(rbind, res)
  }

 result_list <- mclapply(1:4, apply_scheme, mc.cores= detectCores())
  
  for (i in 1:4){
    list_temp[[i + 1]] <- result_list[[i]]
    list_temp[[i + 1]] <- list_temp[[i + 1]][,c(1,7:9,11:16)]
    }
  
  return(list_temp )
}
