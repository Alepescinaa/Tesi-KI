################################################
# Fitting a model with msm over panel data
################################################

library(msm)
library(elect)
library(flexsurv)
library(nhm)

setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/Simulated_data_MM")
load("simulation500_MM_all.RData")
#setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/Simulated_data_MM")
#load("simulation5K_SM_all.RData")

setwd("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/FittingPanelData")
source("hazards_mine.R")


data <- dataset_all_MM_500[[5]] 
n_pats <- nrow(data[[1]])


temp <-data[[2]]
temp$state <- 1

update_state <- function(df) {
  df <- df %>%
    group_by(patient_id) %>%
    mutate(
      state = ifelse(onset == 1, 2, state),
      state = ifelse(row_number() == n() & dead == 1, 3, state),
      #state = ifelse(row_number() == n() & dead == 0, 99, state)
    ) %>%
    ungroup() 
  
  return(df)
}

temp <- update_state(temp)

# =================
# fitting a model
# =================

statetable.msm(state, patient_id, data=temp) #frequency of observation 

Q <- rbind(c(0, 0.2, 0.4), 
           c(0, 0, 0.5),
           c(0, 0, 0)) 


model.msm <- msm(state ~ age,
                 subject = patient_id,
                 data = temp, 
                 qmatrix = Q,
                 covariates = ~ cov1 + cov2 + cov3 ,
                 gen.inits = TRUE,
                 control = list(fnscale = 1000),
                 #censor = 99,
                 deathexact = TRUE)

hazard.msm(model.msm)

# what we find in msm object ?
# estimates <- transition intensities in the log scale
# estimates.t <- transition intensities in natural scale
# qmodel$qmatrix <- intial value of transition intensities 
# Qmatrices <- contains different matrices, component baseline corresponds to estimates.t without cov effects

qmatrix.msm(model.msm) # estimated intenisity matrix on natural scale, as Qmatrices$baseline with ci

# ============
# prediction
# ============
sojourn.msm(model.msm) # mean sojourn time -1/q_ii (before transitioning)

AIC(model.msm)

pmatrix.msm(model.msm, t=10) # P=exp(Q*t) transiton probability matrix
pmatrix.msm(model.msm, t=10, ci="normal", B=100) #obtaining ci of estimates
pmatrix.msm(model.msm, t=10, covariates = list(cov2=0)) # to account for the effect of risk factors over trans prob
pmatrix.msm(model.msm, t=10, covariates = list(cov2=1))

totlos.msm(model.msm, t=30) # total length of stay in each state 
#adds up all periods that a person is predicted to spend in a state, both now and in the future

efpt.msm(model.msm, tostate=3) # expected hitting time to death, well defined we just have one absorbing state

ppass.msm(model.msm, tot=10) # for a person in state r now, this is the probability that before time t, they will
#visit state s at least once

# ==========================
# adding time dependent covs
# ==========================
# let's use elect library to fit a gompertz model by using age as covariate in msm package

# it is better to shift the time dependent covariate
min_age <- min(temp$age)
temp$age <- temp$age - min_age

# let's add baseline indicator
temp <- temp %>%
  group_by(patient_id) %>%
  mutate(bsline = ifelse(row_number() == 1, 1, 0)) %>%
  ungroup()

baseline_data <- temp[temp$bsline==1,]
range(baseline_data$age)

model_age.msm <- msm(state ~ age,
                     subject = patient_id,
                     data = temp, 
                     qmatrix = Q,
                     covariates = ~ cov1 + cov2 + cov3 + age , #elect wants factors as dummy
                     gen.inits=TRUE,
                     control = list(fnscale = 1000),
                     #censor = 99,
                     center= FALSE,# requested by elect
                     deathexact = TRUE) 

# hazard ratios
hazard.msm(model_age.msm) 

# Graph with transition-specific hazard functions derived from an age-dependent model fitted using msm
# change cov values to see how they affect the baseline hazard functions per transition
hazards(model_age.msm, b.covariates = list(age = 0, cov1 = 0, cov2 = 0, cov3 = 0), no.years = 40, age.shift = -min_age)

# The specification of age.max should be such that the probability to survive beyond age.max is assumed to be negligible
# I set all cov effect to zero in order to obtain baseline hazard estimation
elect_model <- elect( x = model_age.msm, b.covariates = list( age = 0, cov1 = 0, cov2 = 0, cov3 = 0),
                      statedistdata = baseline_data, h = 0.1, age.max = 50)

# Life expectancies derived from MLE of model parameters
# defined as integral over time of the transition probability p_ij
elect_model$pnt




#################
# Diagnostic
#################

load("ground_truthMM.RData")

# For each transition defined in trans, the function computes the hazard for that transition
# for each age in the range [age, age + no.years].
haz <- hazards_mine(model_age.msm, b.covariates = list(age = 0, cov1 = 0, cov2 = 0, cov3 = 0), no.years = 40, age.shift = -min_age)

# Assuming this hazards come from fitting a gompertz model I wanna retrieve for each transition 
# shape and rate value, parameters of the distribution
# h(t)=rate*exp(shape*t)
# log(h(t))= log(rate) + shape*t 
# can be seen as y(t)= a+b*t
shape <- numeric()
rate <- numeric()

for (i in 1:3){
  age_grid <- seq(min(temp$age), max(temp$age), length.out = length(unlist(haz[[1]])))
  y <- log(as.numeric(unlist(haz[[i]])))
  reg_model <- lm(y ~ age_grid)
  rate[[i]] <- reg_model$coefficients[1]
  shape[[i]]<- reg_model$coefficients[2] 
}
shape
rate

################################################
# Fitting a model with nhm over panel data
################################################

temp <- data[[2]]
temp$state <- 1

update_state <- function(df) {
  df <- df %>%
    group_by(patient_id) %>%
    mutate(
      state = ifelse(onset == 1, 2, state),
      state = ifelse(row_number() == n() & dead == 1, 3, state)
    ) %>%
    ungroup() 
}
  
temp <- update_state(temp)

tmat_1 <- rbind(c(0,1,2),c(0,0,3),rep(0,3))
tmat_2 <- rbind(c(0,4,5),c(0,0,6),rep(0,3))
tmat_3 <- rbind(c(0,7,8),c(0,0,9),rep(0,3))


temp$patient_id <- as.factor(temp$patient_id) 
temp=as.data.frame(temp)

#here I restrained all covs to have same effect over the same transition

object_nhm <- model.nhm(state ~ age,
                        subject = patient_id,
                        data = temp, 
                        trans = tmat_1,
                        nonh = tmat_1,
                        type = "gompertz",
                        covariates = c("cov1", "cov2", "cov3"),
                        covm = list(cov1= tmat_1, cov2=tmat_2, cov3=tmat_3),
                        #censor.states = c(3),
                        #censor = 99,
                        death = T, 
                        death.states = c(3))


model_nhm <- nhm(object_nhm, 
                 gen_inits = TRUE, 
                 score_test = FALSE, 
                 control = nhm.control(ncores = 4, obsinfo = FALSE))
# we get evidence of nonhomogeneity checking with score test for all transitions
#save(model_nhm, file="model_nhm.RData")


load("model_nhm.RData")
plot(model_nhm, xlim = c(60,100))
plot(model_nhm,what="intensities")


params_nhm <- matrix(0,3,5)
colnames(params_nhm) <- c("rate", "shape", "cov1", "cov2", "cov3")
rownames(params_nhm) <- c("1->2", "1->3", "2->3")
params_nhm[,1] <- model_nhm$par[1:3]
params_nhm[,2] <- model_nhm$par[4:6]
params_nhm[,3] <- model_nhm$par[7:9]
params_nhm[,4] <- model_nhm$par[10:12]
params_nhm[,5] <- model_nhm$par[13:15]


print(model_nhm)

