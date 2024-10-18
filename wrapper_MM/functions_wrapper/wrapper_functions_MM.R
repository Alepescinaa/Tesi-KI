wrapper_functions_MM <- function(data,n_pats){
  
  comp_time <- numeric()
  
  if (scheme==2){
    model_dir <- paste0("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/wrapper_MM/saved_models_scheme2/seed_", seed)
    dir.create(model_dir, showWarnings = FALSE)
  } else if (scheme==3){
    model_dir <- paste0("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/wrapper_MM/saved_models_scheme3/seed_", seed)
    dir.create(model_dir, showWarnings = FALSE)
  } else if (scheme==4){
    model_dir <- paste0("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/wrapper_MM/saved_models_scheme4/seed_", seed)
    dir.create(model_dir, showWarnings = FALSE)
  } else if (scheme==5){
    model_dir <- paste0("/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/wrapper_MM/saved_models_scheme5/seed_", seed)
    dir.create(model_dir, showWarnings = FALSE)
  }
  
  #####################
  # coxph model
  #####################
  
  temp <-  prepare_coxph_flex(data, n_pats)
  model_cox <- vector(mode = "list", length = 3)
  
  time_cox<- system.time({
    for (i in 1:3) {
      model_cox[[i]] <- coxph(Surv(Tstart,Tstop,status) ~ cov1 + cov2 + cov3, data = subset(temp, trans == i))}
  })[3]
  
  comp_time[1] <- as.numeric(round(time_cox,3))
  
  save(model_cox, file = file.path(model_dir,"cox_model.RData"))
  
  ######################
  # flexsurv model
  ######################
  
  temp <- prepare_coxph_flex(data, n_pats)
  fits_gompertz <- vector(mode = "list", length = 3)
  
  time_gomp<- system.time({
    for (i in 1:3) {
      fits_gompertz[[i]] <- flexsurvreg(Surv(Tstart, Tstop, status) ~ cov1 + cov2 + cov3, 
                                        data = subset(temp, trans == i), 
                                        dist = "gompertz")}
  })[3]
  
  comp_time[2] <- as.numeric(round(time_gomp,3))
  
  save(fits_gompertz, file = file.path(model_dir, "flexsurv_model.RData"))
  
  ######################
  # msm model
  ######################
  
  temp <- prepare_msm(data)

  time_msm<- system.time({
    model.msm <- msm(state ~ age,
                     subject = patient_id,
                     data = temp,
                     qmatrix = Q,
                     covariates = ~ cov1 + cov2 + cov3 ,
                     gen.inits = TRUE,
                     control = list(fnscale = 1000, maxit = 500),
                     deathexact = TRUE)
  })[3]

  comp_time[3] <- as.numeric(round(time_msm,3))

  initial_guess_msmage <- qmatrix.msm(model.msm)$estimates
  msm_estimates <- model.msm$estimates.t

  save(model.msm, file = file.path(model_dir, "msm_model.RData"))
  
  #####################
  # msm + age model
  #####################
  
  temp <- prepare_msm(data)

  Q <- rbind(c(0, 1, 1),
             c(0, 0, 1),
             c(0, 0, 0))

  time_msm_age<- system.time({
    model.msm_age <- msm(state ~ age,
                         subject = patient_id,
                         data = temp,
                         qmatrix = initial_guess_age,
                         covariates = ~ cov1 + cov2 + cov3 + age ,
                         gen.inits= TRUE,
                         control = list(fnscale = 1000, maxit = 500),
                         center= TRUE,
                         deathexact = TRUE)
  })[3]

  comp_time[4] <- as.numeric(round(time_msm_age,3))

  save(model.msm_age, file = file.path(model_dir, "model_msm_age.RData"))
  
  ######################
  # nhm model
  ######################
  
  # temp <- prepare_msm(data)
  # 
  # tmat_1 <- rbind(c(0,1,2),c(0,0,3),rep(0,3))
  # tmat_2 <- rbind(c(0,4,5),c(0,0,6),rep(0,3))
  # tmat_3 <- rbind(c(0,7,8),c(0,0,9),rep(0,3))
  # 
  # temp$patient_id <- as.factor(temp$patient_id)
  # temp=as.data.frame(temp)
  # 
  # initial_guess <-  append(msm_estimates, rep(0.5, 3), after = 3)
  # # we have estimates for rate and covs, so we add initial estimate to 0.5 of shape -> not helping
  # 
  # time_nhm <- system.time({
  #   object_nhm <- model.nhm(state ~ age,
  #                           subject = patient_id,
  #                           data = temp,
  #                           trans = tmat_1,
  #                           nonh = tmat_1,
  #                           type = "gompertz",
  #                           covariates = c("cov1", "cov2", "cov3"),
  #                           covm = list(cov1= tmat_1, cov2=tmat_2, cov3=tmat_3),
  #                           death = T,
  #                           death.states = c(3))
  # 
  # 
  #   model_nhm <- nhm(object_nhm,
  #                    gen_inits = FALSE,
  #                    initial = initial_guess,
  #                    score_test = FALSE,
  #                    control = nhm.control(ncores = 1, obsinfo = FALSE, coarsen = T, coarsen.vars = c(1), coarsen.lv = 10))
  # })[3]
  # 
  # #comp.time 500 pats, no coarsening 30 min
  # #comp.time 500 pats, coarsening to 10 values 5min
  # 
  # comp_time[5] <- as.numeric(round(time_nhm,3))
  # 
  # save(model_nhm, file = file.path(model_dir, "model_nhm.RData"))

  ####################
  # imputation model
  ####################
  
  temp <- prepare_imputation(data)
  m <- 30
  type <- "forward"
  
  time_imp <- system.time({
    results_imp<- run_imputation(temp[[1]], temp[[2]], m, type)
    avg_parameters <- results_imp[[1]]
    all_fits <- results_imp[[2]]
  })[3]
  
  comp_time[6] <- as.numeric(round(time_imp,3))
  
  save(results_imp, file = file.path(model_dir, "results_imp.RData"))
  
  save(comp_time, file = file.path(model_dir, "computational_time.RData"))
  
  
}