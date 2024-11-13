run_imputation <- function(data, data_visits, m, type){
  
  schemeA_data <- data
  schemeA_visits <- data_visits
  
  temp <- schemeA_data  
  temp$midage <- 0
  temp$midage <- (schemeA_data$last_bfo + schemeA_data$onset_age)/2
  temp$onset_age <- ifelse(temp$onset==1,temp$midage,temp$onset_age )
  mod_onset <- coxph(Surv(age,onset_age,onset)~cov1+cov2+cov3, data = temp)
  onsethaz<-basehaz(mod_onset,center=FALSE) 
  onsetpar<-mod_onset$coefficients
  
  temp$start<-temp$age
  temp$start <- ifelse(temp$onset==1,temp$midage, temp$start )
  mod_death<- coxph(Surv(start,death_time,dead)~cov1+cov2+cov3+onset, data = temp)
  deathhaz<-basehaz(mod_death,center=FALSE)
  deathpar<-mod_death$coefficients
  
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
    remove(tzero)
    
    
    combhaz<-merge(onsethaz,deathhaz,by.x="time", by.y="time",all.x=TRUE, all.y=TRUE, suffixes = c(".12",".13"))
    #merging hazards from both models, since they are computed on different transition times we might have NA values
    if (is.na(combhaz$hazard.12[1])) combhaz$hazard.12[1]<-0
    if (is.na(combhaz$hazard.13[1])) combhaz$hazard.13[1]<-0
    
    #avoiding NA
    for (i in 2: nrow(combhaz))
    {
      if (is.na(combhaz$h.12[i]))
      {
        combhaz$h.12[i]<-combhaz$h.12[i-1]
        combhaz$hazard.12[i]<-combhaz$hazard.12[i-1]+combhaz$h.12[i-1]*(combhaz$time[i]-combhaz$time[i-1])
      }
      if (is.na(combhaz$h.13[i]))
      {
        combhaz$h.13[i]<-combhaz$h.13[i-1]
        combhaz$hazard.13[i]<-combhaz$hazard.13[i-1]+combhaz$h.13[i-1]*(combhaz$time[i]-combhaz$time[i-1])
      }
    }
    
    demstatus<-matrix(nrow=nrow(schemeA_data),ncol=m,rep(schemeA_data$onset,m))
    # matrix representing dementia status for each patient for each imputated dataset (set to 1 if dementia observed in original dataset)
    dementage<-matrix(nrow=nrow(schemeA_data),ncol=m,0)
    # matrix representing dementia onset age for each patient for each imputated dataset (in this case is censored in orignal dataset, so set to zero for everyone)
    c1<-matrix(nrow=nrow(schemeA_data),ncol=m,runif(nrow(schemeA_data)*m)) # matrix patients* n_imputed_datasets filled of numbers sampled in [0,1]
    c2<-matrix(nrow=nrow(schemeA_data),ncol=m,runif(nrow(schemeA_data)*m))
    covar<-as.matrix(schemeA_data[,c("cov1","cov2", "cov3")]) #covariates values for each patient
    lp.12<-covar %*% as.vector(t(onsetpar[c("cov1","cov2", "cov3")])) #covariates* theirs estimated coefficients of dementia onset model
    lp.13<-covar %*% as.vector(t(deathpar[c("cov1","cov2", "cov3")]))  #covariates* theirs estimated coefficients of death model
    delta<-deathpar["onset"]
    
    
    for (i in 1:nrow(schemeA_data)){
      time<-combhaz$time
      S1<-exp(-combhaz$hazard.12*exp(lp.12[i])-combhaz$hazard.13*exp(lp.13[i])) #survival function for living dementia free
      S2<-exp(-combhaz$hazard.13*exp(lp.13[i]+delta)) #survival function for living
      p<-S1*combhaz$h.12/S2 # instantaneous probability of transitioning to dementia conditional on survival upon that time
      
      P<-p   
      P[1]<-0 #P is gonna be cumulative probability of developing dementia over time
      for (j in 2:length(P)) P[j]<-min(1.0,P[j-1]+p[j-1]*(time[j]-time[j-1]))
      F<-data.frame(cbind(time,S1,S2,p,P))
      # names(F)<-c("time","S1","S2","p","P")
      remove(time,S1,S2,p,P)
      
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
        a<-patient_data$age[n_visits]
        b<-patient_data$death_time[1]
        d<-patient_data$dead[1]
        F<-F[F$time>a & F$time<b,] #the onset time is restricted to interval (last visit, death/censoring time)
        nF<-nrow(F)
        
        if (nF>0)# if there are observation times in that interval
        {     
          minP<-F$P[1]     
          maxP<-F$P[nF]
          maxt<-F$time[nF]
          S2t<-F$S2[nF] #surival prob of not dying at maximum timepoint in a,b
          S1t<-F$S1[nF] #surival prob of not developing dementia at maximum timepoint in a,b
          S2h23P<-S2t*exp(lp.13[i]+d*delta)*(maxP-minP) #adjust survival probability (extra weight given to disease if patients die) 
          # basically we are tring to give more probability of being extracted as ill patient to those that are dead since
          # dead people could have died because of unobserved dementia
          pD<-S2h23P/(S1t+S2h23P)    # probability of dementia onset during the interval, given survival without dementia up to the final time point  
          
          if (minP<maxP) F$q<-(F$P-minP)/(maxP-minP)  #normalized cum probabilities of developing dementia, needed to set onset-age
          if (minP>=maxP) F$q<-0
        }
        if (nF<=0) pD<-(a<b)*mean(schemeA_data$onset)# if there are observation times in that interval probability of dementia onset 
        #set as the proportion of observed demented people
        
        for (j in 1:m) #creating m imputated datasets
        {
          if (c2[i,j]<=pD) #if  probability of dementia onset is more than a threshold dementia status=1
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
    
    #fitting the parametric model for each imputed dataset and averaging estimates
    for (j in 1:m){
      temp <- schemeA_data
      
      dem_generated <- which(demstatus[,j]==1)
      dem_observed <- as.numeric(temp$patient_id[temp$onset==1])
      dem_to_add <- setdiff(dem_generated, dem_observed)
      
      temp$onset[dem_to_add] <- 1
      temp$onset_age[dem_generated] <- dementage[dem_generated,j]
      
      all_fits[[j]] <- fit_model(temp,type)
      
    }
    
    averaged_params <- averaging_params(all_fits,m)

  return(list(averaged_params = averaged_params, all_fits = all_fits))
}