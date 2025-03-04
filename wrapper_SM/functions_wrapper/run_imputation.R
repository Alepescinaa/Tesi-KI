run_imputation <- function(data, data_visits, m, type, cores_imp){
  
  schemeA_data <- data
  schemeA_visits <- data_visits
  
  # we first impute the onset age for those having dementia observed at the middle point between 
  # previous visit and visit of the diagnosis
  schemeA_data$midage <- 0
  schemeA_data$midage <- (schemeA_data$last_bfo + schemeA_data$onset_age)/2
  schemeA_data$onset_age <- ifelse(schemeA_data$onset==1,schemeA_data$midage,schemeA_data$onset_age )
  schemeA_data <- schemeA_data %>%
    mutate(midage=NULL,last_bfo=NULL)
  schemeA_data$entry <- 0
  schemeA_data$entry[schemeA_data$onset==1] <- schemeA_data$onset_age[schemeA_data$onset==1] 

  tmat <- mstate::transMat(x = list(c(2, 3),c(3),c()), names = c("Dementia-free","Dementia", "Death")) 
  n_trans <- max(tmat, na.rm = TRUE)
  
  data_long <- msprep(data = schemeA_data, trans = tmat, 
                      time = c(NA, "onset_age", "death_time"), 
                      status = c(NA, "onset", "dead"), 
                      keep = c("age","cov1", "cov2", "cov3"),
                      id="patient_id")

  
  data_long$Tstart[data_long$trans<3] <- data_long$Tstart[data_long$trans<3] +data_long$age[data_long$trans<3]
  data_long$time <- data_long$Tstop-data_long$Tstart
  data_long$entry <- 0
  data_long$entry[data_long$trans==3] <- data_long$Tstart[data_long$trans==3]
  data_long <- expand.covs(data_long,c("cov1","cov2", "cov3", "entry"))
  data_long$Tstart2 <-  data_long$Tstart
  data_long$Tstop2 <-  data_long$Tstop
  data_long$Tstart2[data_long$trans==3] <- 0
  data_long$Tstop2[data_long$trans==3] <- data_long$time[data_long$trans==3]
  data_long <-data_long %>%
    mutate(Tstart=Tstart2, Tstop=Tstop2, Tstart2=NULL, Tstop2=NULL)
  model_cox<- coxph(Surv(Tstart, Tstop ,status) ~ cov1.1 + cov2.1 + cov3.1 + cov1.2 + cov2.2 + cov3.2 + cov1.3 + cov2.3 + cov3.3 + entry.1 + entry.2 + entry.3 + strata(trans), data = data_long, method="breslow")
  all_haz <- basehaz(model_cox,center=FALSE) 
  
  onsethaz <- all_haz[all_haz$strata=="trans=1",1:2]
  onsetpar<-model_cox$coefficients[1:3]
  
  deathhaz <- all_haz[all_haz$strata=="trans=2",1:2]
  deathpar<-model_cox$coefficients[4:6]
  
  deathhaz_dem <- all_haz[all_haz$strata=="trans=3",1:2]
  deathpar_dem<-model_cox$coefficients[c(7:9,12)]

  
    # hazards computation
    tzero<-data.frame(matrix(0,1,2))
    names(tzero)<-c("time","hazard")
    onsethaz<-rbind(tzero,onsethaz) 
    onsethaz$h<-onsethaz$hazard
    onsethaz$h[nrow(onsethaz)]<-NA
    for (i in 1:(nrow(onsethaz)-1))
      onsethaz$h[i]<-(onsethaz$hazard[i+1]-onsethaz$hazard[i])/(onsethaz$time[i+1]-onsethaz$time[i])
    
    deathhaz<-rbind(tzero,deathhaz)
    deathhaz$h<-deathhaz$hazard
    deathhaz$h[nrow(deathhaz)]<-NA
    for (i in 1:(nrow(deathhaz)-1))
      deathhaz$h[i]<-(deathhaz$hazard[i+1]-deathhaz$hazard[i])/(deathhaz$time[i+1]-deathhaz$time[i])
    
    
    deathhaz_dem<-rbind(tzero,deathhaz_dem)
    deathhaz_dem$h<-deathhaz_dem$hazard
    deathhaz_dem$h[nrow(deathhaz_dem)]<-NA
    for (i in 1:(nrow(deathhaz_dem)-1))
      deathhaz_dem$h[i]<-(deathhaz_dem$hazard[i+1]-deathhaz_dem$hazard[i])/(deathhaz_dem$time[i+1]-deathhaz_dem$time[i])
    remove(tzero)
   
   
    combhaz<-merge(onsethaz,deathhaz,by.x="time", by.y="time",all.x=TRUE, all.y=TRUE, suffixes = c(".01",".02"))
   
    #merging hazards from both models, since they are computed on different transition times we might have NA values
    if (is.na(combhaz$hazard.01[1])) combhaz$hazard.01[1]<-0
    if (is.na(combhaz$hazard.02[1])) combhaz$hazard.02[1]<-0
    #if (is.na(combhaz$hazard.02[1])) combhaz$hazard.12[1]<-0
    
    #avoiding NA
    for (i in 2: nrow(combhaz))
    {
      if (is.na(combhaz$h.01[i]))
      {
        combhaz$h.01[i]<-combhaz$h.01[i-1]
        combhaz$hazard.01[i]<-combhaz$hazard.01[i-1]+combhaz$h.01[i-1]*(combhaz$time[i]-combhaz$time[i-1])
      }
      if (is.na(combhaz$h.02[i]))
      {
        combhaz$h.02[i]<-combhaz$h.02[i-1]
        combhaz$hazard.02[i]<-combhaz$hazard.02[i-1]+combhaz$h.02[i-1]*(combhaz$time[i]-combhaz$time[i-1])
      }
    }
  
    demstatus<-matrix(nrow=nrow(schemeA_data),ncol=m,rep(schemeA_data$onset,m))
    # matrix representing dementia status for each patient for each imputated dataset (set to 1 if dementia observed in original dataset)
    dementage<-matrix(nrow=nrow(schemeA_data),ncol=m,0)
    # matrix representing dementia onset age for each patient for each imputated dataset (in this case is censored in orignal dataset, so set to zero for everyone)
    c1<-matrix(nrow=nrow(schemeA_data),ncol=m,runif(nrow(schemeA_data)*m)) # matrix patients* n_imputed_datasets filled of numbers sampled in [0,1]
    c2<-matrix(nrow=nrow(schemeA_data),ncol=m,runif(nrow(schemeA_data)*m))
    covar<-as.matrix(schemeA_data[,c("cov1","cov2", "cov3")]) #covariates values for each patient
    lp.01<-covar %*% as.vector(t(onsetpar[c("cov1.1","cov2.1", "cov3.1")])) #covariates* theirs estimated coefficients of dementia onset model
    lp.02<-covar %*% as.vector(t(deathpar[c("cov1.2","cov2.2", "cov3.2")]))  #covariates* theirs estimated coefficients of death model
    lp.12<-covar %*% as.vector(t(deathpar_dem[c("cov1.3","cov2.3", "cov3.3")]))  #covariates* theirs estimated coefficients of death model having dementia
    delta<-deathpar_dem["entry.3"]
    
    
    for (i in 1:nrow(schemeA_data)){

      time<-combhaz$time
      S1<-exp(-combhaz$hazard.01*exp(lp.01[i])-combhaz$hazard.02*exp(lp.02[i])) #survival function for living dementia free
      S2<-exp(-combhaz$hazard.02*exp(lp.02[i])) #survival function for living (conditional you are healthy)
      S2_dem<-exp(-deathhaz_dem$hazard*exp(lp.12[i])) #survival function for living (conditional you are demented)
      
      p<-S1*combhaz$h.01/S2 # instantaneous probability of transitioning to dementia conditional on survival upon that time
      
      P<-p   
      P[1]<-0 #P is gonna be cumulative probability of developing dementia over time
      for (j in 2:length(P)) P[j]<-min(1.0,P[j-1]+p[j-1]*(time[j]-time[j-1]))
      F<-data.frame(cbind(time,S1,S2,p,P))
      dem_df <- data.frame("time"=deathhaz_dem$time, S2_dem)
    
  
      remove(time,S1,S2, S2_dem,p,P)
      
      patient_data <- schemeA_visits[schemeA_visits$patient_id==i,]
      
      #generate onset time for those having dementia observed
      
      if (any(patient_data$onset==1)){
        index <-which(patient_data$onset==1)[1] #first visit observing dementia
        a <- patient_data$age[index-1]
        b <- patient_data$age[index]
        F<-F[F$time>a & F$time<b,]
        nF<-nrow(F) 
        
        if (nF>0) # if there are observation times in that interval
        {    
          minP<-F$P[1]
          maxP<-F$P[nF] 
          if (minP<maxP) #if the probability of transitioning increases over the interval
          {
            F$q<-(F$P-minP)/(maxP-minP)  #normalized probabilities
            for (j in 1:m) dementage[i,j]<-max(F[F$q<=c1[i,j],]$time) #onset is set at maximum time for which probabiltiy is less of
            #a chose treshold sampled from a uniform
          }
          if (minP>=maxP)
            dementage[i,]<-(a+b)/2
        }
        if (nF<=0) #if no observation times are available we set onset in the middle point
          dementage[i,]<-(a+b)/2
      }
      
      #generate onset value and onset time for censored patients
      if (!any(patient_data$onset==1)) {
        n_visits <- length(patient_data$visits)
        a<-patient_data$age[n_visits] # age last visit
        b<-patient_data$death_time[1] # age lost/death
        d<-patient_data$dead[1] #indicator od death 
        F<-F[F$time>a & F$time<b,] #the onset time is restricted to interval (last visit, death/censoring time)
        int_length <- b-a
        # in dem_df I have the value of the survival function for living associated to an amount of time spent
        # in dementia status, I am gonna consider those values closer to the time spent from
        # last visit to death for the cosnidered patient, ie. the time he spent in dementia status
        # assuming he had developed dementia after last visit
        dem_df <- dem_df[order(abs(dem_df$time - int_length)), ][1, ]
        nF<-nrow(F)
        
        if (nF>0)# if there are observation times in that interval
        {     
          minP<-F$P[1]     
          maxP<-F$P[nF] 
          maxt<-F$time[nF]
          S2t<-F$S2[nF] #survival fun of not dying at maximum timepoint in a,b
          S1t<-F$S1[nF] #survival fun of not developing dementia at maximum timepoint in a,b
          S2_dem_a <- dem_df$S2_dem #survival fun of not dying having had dementia for time = b-a
          S2h23P<-(S2t*exp(lp.02[i])*(maxP-minP) +d*(1-S2_dem_a)*exp(lp.12[i]+delta*a))*(maxP-minP) #adjusted probability of developing dem (extra weight given to disease if patients die) 
          # basically we are trying to give more probability of being extracted as ill patient to those that are dead since
          # dead people could have died because of unobserved dementia, so we are adding the survival function 
          # associated to the time spent in dementia status adjusted for the entry in dementia status
          pD<-S2h23P/(S1t+S2h23P)   # probability of dementia onset during the interval, given survival without dementia up to the final time point  
          
          
          if (minP<maxP) F$q<-(F$P-minP)/(maxP-minP)  #normalized cum probabilities of developing dementia, needed to set onset-age
          if (minP>=maxP) F$q<-0
        }
        if (nF<=0) pD<-(a<b)*mean(schemeA_data$onset)# if there are no observation times in that interval probability of dementia onset 
        #set as the proportion of observed demented people
        
        for (j in 1:m) #creating m imputated datasets
        {
          if (c2[i,j]<=pD) #if probability of dementia onset is more than a threshold dementia status=1
            # and onset age computed as before: at maximum time for which cum probability of dementia is less of
            #a chosen threshold sampled from a uniform or in middle point
          {
            demstatus[i,j]<-1 
            if (nF>0) dementage[i,j]<-max(F[F$q<=c1[i,j],]$time)
            if (nF<=0) dementage[i,j]<-(a+b)/2
          }
          if (c2[i,j]>pD) #if  probability of dementia onset is less than a threshold dementia status=0
          {
            demstatus[i,j]<-0
            dementage[i,j]<-patient_data$onset_age[1] # we set onset_age at the censoring value
            
            
            
          }
        }
      }
    }
    
  
    
    all_fits <- vector(mode = "list", length = m)
    
    plan(multisession, workers = cores_imp)
    
    # Run loop in parallel
    all_fits <- future_lapply(1:m, function(j) {
      temp <- schemeA_data
      
      dem_generated <- which(demstatus[,j] == 1)
      dem_observed <- as.numeric(temp$patient_id[temp$onset == 1])
      dem_to_add <- setdiff(dem_generated, dem_observed)
      
      temp$onset[dem_to_add] <- 1
      temp$onset_age[dem_generated] <- dementage[dem_generated, j]
      
      fit_model(temp, type)  
    })
    
    # Reset plan to sequential after parallel execution
    plan(sequential)
    
    averaged_params <- averaging_params(all_fits,m)

  return(list(averaged_params = averaged_params, all_fits = all_fits))
}

