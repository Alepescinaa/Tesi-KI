################################################
# Fitting a coxph MM model over Panel Data
################################################

library (survival)
library(survminer)
setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/Simulated_data_MM")
load("simulation5K_MM_all.RData")

setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/FittingPanelData")


data <- dataset_all_MM_2K[[2]] 
n_pats <- nrow(data[[1]])

scheme_data <- data[[2]]

original_index <- unique(scheme_data$patient_id)
scheme_data$last_bfo <- numeric()
for (i in original_index){
  index <- which(scheme_data$onset[scheme_data$patient_id==i]==1)[1]
  scheme_data$last_bfo[scheme_data$patient_id==i] <- scheme_data$age[scheme_data$patient_id==i][index-1]
}
for (i in original_index){
  if(any(scheme_data$onset[scheme_data$patient_id==i]==1))
    scheme_data$onset[scheme_data$patient_id==i] <- 1
}

scheme_data <- scheme_data[order(scheme_data$patient_id),]

row_id <- scheme_data %>%
  group_by(patient_id) %>%  
  summarise(nrows = n())

scheme_data$patient_id <- rep(1:n_pats, times=as.numeric(row_id$nrows))

scheme_data <- scheme_data %>%
  group_by(patient_id) %>%
  slice(1) %>%
  ungroup()

tmat <- mstate::transMat(x = list(c(2, 3),c(3),c()), names = c("Dementia-free","Dementia", "Death")) 
n_trans <- max(tmat, na.rm = TRUE)

data_long <- msprep(data = scheme_data, trans = tmat, 
                       time = c(NA, "onset_age", "death_time"), 
                       status = c(NA, "onset", "dead"), 
                       keep = c("age","cov1", "cov2", "cov3"),
                       id="patient_id")

data_long$Tstart[data_long$trans<3] <- data_long$Tstart[data_long$trans<3] +data_long$age[data_long$trans<3]
data_long$time <- data_long$Tstop-data_long$Tstart

#============================
# coxph model per transition
#============================
model_cox <- vector(mode = "list", length = n_trans)

for (i in 1:3) {
  model_cox[[i]] <- coxph(Surv(Tstart,Tstop,status) ~ cov1 + cov2 + cov3, data = subset(data_long, trans == i))}
 
model_cox[[1]]
model_cox[[2]]
model_cox[[3]]

par(mfrow=c(3,1))

surv_fit <- survfit(model_cox[[1]],newdata = subset(data_long, trans == 1))

plot(surv_fit, 
     xlim = c(60, 100),
     ylim= c(0, 1),
     xlab = "Age", 
     ylab = "Survival Probability", 
     main = "Survival Curves by Dementia v.1", 
     col = c("blue"), 
     lty = 1:2,                
     lwd = 1)  

#====================================
# coxph coherent with imputation code
#====================================

# imputation of onset age at middle age and considering just right censoring 
temp <- scheme_data  
temp$midage <- 0
temp$midage <- (scheme_data$last_bfo + scheme_data$onset_age)/2
temp$onset_age <- ifelse(temp$onset==1,temp$midage,temp$onset_age )
mod_onset_rc <- coxph(Surv(age,onset_age,onset) ~cov1+cov2+cov3, data = temp)
mod_onset_rc
AIC(mod_onset_rc)

# VERSION TO CONSIDER INTERVAL CENSORING
# temp <- scheme_data
# temp$onset <- ifelse(temp$onset==1, 3, temp$onset) #to specify that onset age is interval censored ????
# mod_onset_ic <- coxph(Surv(age,onset_age,onset, type=c("interval2"))~cov1+cov2+cov3, data = temp)
# mod_onset_ic
# AIC(mod_onset_ic) 

#===============
# death model
#===============

# coherent with first onset model
temp <- scheme_data  
temp$midage <- 0
temp$midage <- (scheme_data$last_bfo + scheme_data$onset_age)/2
temp$start<-temp$age
temp$start <- ifelse(temp$onset==1,temp$midage, temp$start )
mod_death_rc<- coxph(Surv(start,death_time,dead)~cov1+cov2+cov3+strata(onset), data = temp)
mod_death_rc
AIC(mod_death_rc)



# ===========
# plot curves
# ===========
par(mfrow=c(2,1))

# surv_fit_onset <- survfit(mod_onset_rc)
# ggsurv_onset_rc <- ggsurvplot(surv_fit_onset, data = temp,
#                               conf.int = TRUE,        
#                               pval = TRUE,            
#                               risk.table = TRUE,      
#                               ggtheme = theme_minimal())  

surv_fit <- survfit(mod_onset_rc0)
plot(surv_fit, 
     xlim = c(60, 100),
     ylim= c(0, 1),
     xlab = "Age", 
     ylab = "Survival Probability", 
     main = "Survival Curves by Dementia v.1", 
     col = c("blue"), 
     lty = 1:2,                
     lwd = 1)    

surv_fit <- survfit(mod_death_rc)
plot(surv_fit, 
     xlim = c(60, 100), 
     xlab = "Age", 
     ylab = "Survival Probability", 
     main = "Survival Curves by Death v.2", 
     col = c("blue", "red"), 
     lty = 1:2,                
     lwd = 2)                 
