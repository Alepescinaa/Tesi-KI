apply_right_censoring <- function(temp){
  right_censoring <- runif (n=n_pats, min=0, max=20) 
  rc_id <- data_frame(rc = right_censoring,
                      id = 1:n_pats)
  temp$right_censoring[temp$final == 1] <- rc_id$rc[match(temp$patient_id[temp$final == 1], rc_id$id)] +   temp$time_start[temp$final == 1]
  if (all(temp$right_censoring[temp$final == 1]>temp$time_start[temp$final == 1])==FALSE)
    print(error)
  temp$dead[temp$final==1] <- if_else(temp$time_stop[temp$final==1] > temp$right_censoring[temp$final==1],0,1)
  temp$death_time[temp$final==1]<-if_else(temp$time_stop[temp$final==1]<temp$right_censoring[temp$final==1],temp$time_stop[temp$final==1],temp$right_censoring[temp$final==1]) 
  return(temp)
}