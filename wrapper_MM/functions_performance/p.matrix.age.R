p.matrix.age <- function(x = NULL, age1 = 1, age2 = 0, int=0.1, covariates = "mean", ci = c("none", 
                                                                                            "normal", "bootstrap"), cl = 0.95, B = 1000, cores = NULL, 
                         qmatrix = NULL, ...){
  
  age <- seq(age1,age2,by=int)
  d <- dim(x$qmodel$imatrix)
  prob <- array(NA,dim=c(d[1],d[2],length(age)))
  if(ci!="none"){
    L <- array(NA,dim=c(d[1],d[2],length(age)))
    U <- array(NA,dim=c(d[1],d[2],length(age)))
  }
  covs <- list(do.call(c, list(list(age=age1), covariates)))
  for(a in 2:length(age)){
    cov_age <- do.call(c, list(list(age=age[a]), covariates))
    covs[[a]] <- cov_age
  }
  
  
  
  for (a in 2:length(age)){
    start <- age1+int
    times <- seq(start,age[a],by=int)
    
    prob_obj <- pmatrix.piecewise.msm(x,
                                      t1 = min(age),
                                      t2=age[a],
                                      times = times,
                                      covariates =covs[1:(length(times)+1)] , 
                                      ci, cl ,B, cores , 
                                      qmatrix , ...)
    if(ci=="none"){
      prob[,,a-1] <- prob_obj
    }else{
      prob[,,a-1] <- prob_obj$estimate
      L[,,a-1] <- prob_obj$L
      U[,,a-1] <- prob_obj$U
    }
    
  }
  
  if(ci=="none"){
    prob_molt <- data.frame(trans=rep(paste0(rep(paste0("from",1:d[1]),d[1]),rep(paste0("to",1:d[2]),each=d[2])),length(age)),
                            age=rep(age,each=d[1]*d[2]),
                            prob=as.numeric(prob))
  }else{
    prob_molt <- data.frame(trans=rep(paste0(rep(paste0("from",1:d[1]),d[1]),rep(paste0("to",1:d[2]),each=d[2])),length(age)),
                            age=rep(age,each=d[1]*d[2]),
                            prob=as.numeric(prob),
                            upper=as.numeric(U),
                            lower=as.numeric(L))
  }
  
  return(prob_molt)
}