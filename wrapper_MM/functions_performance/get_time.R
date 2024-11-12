get_time <- function(prob,ci="none",int=0.1){
  if(ci=="none"){
    prob %<>% group_by(trans)%>%  summarise(time=sum(prob*int,na.rm=T))
    
  }else{
    prob %<>% group_by(trans)%>%  summarise(time=sum(prob*int,na.rm=T),
                                            time_lower=sum(lower*int,na.rm=T),
                                            time_upper=sum(upper*int,na.rm=T))
  }
  return(prob)
}