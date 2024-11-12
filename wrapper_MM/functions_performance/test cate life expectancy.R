library(msm)
library(tidyverse)
library(magrittr)
library(haven)
library(msmtools)
library(splines)
source("fun prob msm.R")
dat <- read_dta("data\\dataset_II_dementia_CIND.dta")

trans_q <-
  rbind(c(0, 0.1, 0, 0.05),
        c(0.05, 0, 0.1, 0.06),
        c(0, 0, 0, 0.30),
        c(0, 0, 0, 0))

qcrude <-
  crudeinits.msm(statecens_ ~ age,
                 subject = lopnr,
                 data = dat,
                 qmatrix = trans_q,
                 censor = 5,
                 censor.states = c(1,2))

msm_6_ptau_dic <- msm(
  statecens_ ~ age,
  subject = lopnr,
  data=dat,
  qmatrix = qcrude,
  censor = 5,
  gen.inits = T,
  censor.states = c(1,2),
  covariates = list("1-2" = ~ age + sex + educ_2 + educ_3 + ptau_high,  
                    "2-1" = ~ age + sex + educ_2 + educ_3 + ptau_high,  
                    "2-3" = ~ age + sex + educ_2 + educ_3 + ptau_high,
                    "1-4" = ~ age,
                    "2-4" = ~ age,
                    "3-4" = ~ age),
  exacttimes = F,
  deathexact = T,
  death = 4,
  method = "BFGS",
  control = list(fnscale = 9000, maxit = 15000))
  
  

test <- p.matrix.age(msm_6_ptau_dic,
             covariates = list(
                               sex=0.6,
                               ptau_high=0,
                               educ_2=0,
                               educ_3=0),
             age1 = 60,
             age2 = 100,
             int=0.1,
             ci="none",
             cores=4)


ggplot(test)+
  #geom_ribbon(aes(age,ymin=lower,upper=upper),alpha=0.3)
  geom_line(aes(age,prob))+
  facet_wrap(~trans,scales = "free_y")


  
get_time(test,ci="none")

test <- p.matrix.age(msm_6_ptau_dic,
                     covariates = list(
                       sex=0.6,
                       ptau_high=0,
                       educ_2=0,
                       educ_3=0),
                     age1 = 60,
                     age2 = 100,
                     int=2,
                     ci="normal",
                     cores=4)

ggplot(test)+
  geom_ribbon(aes(age,ymin=lower,upper=upper),alpha=0.3)
  geom_line(aes(age,prob))+
  facet_wrap(~trans,scales = "free_y")

  
get_time(test,ci="normal",int=2)