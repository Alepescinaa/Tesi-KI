wrapper_functions_SM <- function(data,n_pats,seed,cores_nhm){
  
  comp_time <- numeric()
  setwd(here())

  if (n_pats==500){
    if (scheme==2){
      model_dir <- paste0("wrapper_SM/results_500/saved_models_scheme2/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==3){
      model_dir <- paste0("wrapper_SM/results_500/saved_models_scheme3/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==4){
      model_dir <- paste0("wrapper_SM/results_500/saved_models_scheme4/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  
    } else if (scheme==5){
      model_dir <- paste0("wrapper_SM/results_500/saved_models_scheme5/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  }
  } else if (n_pats==2000){
    if (scheme==2){
      model_dir <- paste0("wrapper_SM/results_2K/saved_models_scheme2/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==3){
      model_dir <- paste0("wrapper_SM/results_2K/saved_models_scheme3/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==4){
      model_dir <- paste0("wrapper_SM/results_2K/saved_models_scheme4/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  
    } else if (scheme==5){
      model_dir <- paste0("wrapper_SM/results_2K/saved_models_scheme5/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  }
  } else if (n_pats==5000){
    if (scheme==2){
      model_dir <- paste0("wrapper_SM/results_5K/saved_models_scheme2/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==3){
      model_dir <- paste0("wrapper_SM/results_5K/saved_models_scheme3/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==4){
      model_dir <- paste0("wrapper_SM/results_5K/saved_models_scheme4/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  
    } else if (scheme==5){
      model_dir <- paste0("wrapper_SM/results_5K/saved_models_scheme5/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  }
  } else if (n_pats==10000){
    if (scheme==2){
      model_dir <- paste0("wrapper_SM/results_10K/saved_models_scheme2/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==3){
      model_dir <- paste0("wrapper_SM/results_10K/saved_models_scheme3/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==4){
      model_dir <- paste0("wrapper_SM/results_10K/saved_models_scheme4/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  
    } else if (scheme==5){
      model_dir < paste0("wrapper_SM/results_10K/saved_models_scheme5/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  }
  }

  n_pats <- length(unique(data$patient_id))
  
  #####################
  # coxph model (ignores ic, accounts sm)
  #####################
  
  temp <-  prepare_coxph(data, n_pats)
  
  error <- F
  
  time_coxph <- system.time({
    tryCatch({
      model_cox<- coxph(Surv(Tstart, Tstop ,status) ~ cov1.1 + cov2.1 + cov3.1 + cov1.2 + cov2.2 + cov3.2 + cov1.3 + cov2.3 + cov3.3 + strata(trans), data = temp, method="breslow")
    },
    error = function(e) {
      print(paste("Error during model fitting:", e$message))
      error <<- TRUE
    })
  })[3]
  
  if (error) {
    print(paste("No cox convergence for seed:", seed))
    model_cox <- NULL
  } else {
    print("Model fitted successfully.")
  }
  
  
  comp_time[1] <- as.numeric(round(time_coxph,3))
  
  if (!is.null(model_cox)) {
    save(model_cox, file = file.path(model_dir, "cox_model.RData"))
  } else {
    print("Model Cox is NULL; not saving.")
  }
  
  gc()
  
 
  
  ######################
  # flexsurv model (ignores ic, accounts sm)
  ######################
  
  temp <- prepare_flex(data, n_pats)
  fits_gompertz <- vector(mode = "list", length = 3)

  error <- F
  
  time_gomp <- system.time({
    tryCatch({
      for (i in 1:2){
        fits_gompertz[[i]] <- flexsurvreg(Surv(Tstart, Tstop, status) ~ cov1 + cov2 + cov3,
                                        data = subset(temp, trans == i),
                                        dist = "gompertz")
        }
      fits_gompertz[[3]] <- flexsurvreg(Surv(time, status) ~ cov1 + cov2 + cov3,  
                                        data = subset(temp, trans == 3),
                                        dist = "gompertz")
    },
    error = function(e) {
      print(paste("Error during model fitting:", e$message))
      error <<- TRUE
    })
  })[3]
  
  if (error) {
    print(paste("No gomp convergence for seed:", seed))
    fits_gompertz <- NULL
  } else {
    print("Model fitted successfully.")
  }
  
  
  comp_time[2] <- as.numeric(round(time_gomp,3))
  
  if (!is.null(fits_gompertz)) {
    save(fits_gompertz, file = file.path(model_dir, "flexsurv_model.RData"))
  } else {
    print("Model gomp is NULL; not saving.")
  }
  
  gc()
  
  
  # for nhm parameters are in the following order rate1,rate2,rate3, shape1, shape2, shape3, covs...
  # initial_guess <- c( fits_gompertz[[1]]$coefficients[2], fits_gompertz[[2]]$coefficients[2],
  #                     fits_gompertz[[3]]$coefficients[2], fits_gompertz[[1]]$coefficients[1],
  #                     fits_gompertz[[2]]$coefficients[1], fits_gompertz[[3]]$coefficients[1],
  #                     0,0,0,0,0,0,0,0,0)

  

  ######################
  # nhm model (accounts ic, ignores sm)
  ######################
  
  #would be necessary coarsening more and add more splits wrt markovian case

  # temp <- prepare_msm(data)
  # 
  # tmat_1 <- rbind(c(0,1,2),c(0,0,3),rep(0,3))
  # tmat_2 <- rbind(c(0,4,5),c(0,0,6),rep(0,3))
  # tmat_3 <- rbind(c(0,7,8),c(0,0,9),rep(0,3))
  # 
  # temp$patient_id <- as.factor(temp$patient_id)
  # temp=as.data.frame(temp)
  # 
  # 
  # find_splits <- function(age) {
  #   quantiles <- quantile(age, probs = seq(0, 1, 0.01))
  #   return(quantiles[-c(1, length(quantiles))])
  # }
  # 
  # split_points <- find_splits(temp$age)
  # 
  # 
  # error <- F
  # 
  # time_nhm <- system.time({
  #   tryCatch({
  #     object_nhm <- model.nhm(state ~ age,
  #                             subject = patient_id,
  #                             data = temp,
  #                             trans = tmat_1,
  #                             nonh = tmat_1,
  #                             type = "gompertz",
  #                             covariates = c("cov1", "cov2", "cov3"),
  #                             covm = list(cov1= tmat_1, cov2=tmat_2, cov3=tmat_3),
  #                             censor.states = c(1,2),
  #                             censor = 99,
  #                             death = T,
  #                             death.states = c(3))
  # 
  # 
  #     model_nhm <- nhm(object_nhm,
  #                      gen_inits = T,
  #                      #initial = initial_guess,
  #                      score_test = FALSE,
  #                      control = nhm.control(ncores = cores_nhm, obsinfo = FALSE, coarsen = T, coarsen.vars = c(1), coarsen.lv = 5, splits = split_points, rtol=1e-4, atol=1e-4))
  #   },
  #   error = function(e) {
  #     print(paste("Error during model fitting:", e$message))
  #     error <<- TRUE
  #   })
  # })[3]
  # 
  # if (error) {
  #   print(paste("No nhm convergence for seed:", seed))
  #   model_nhm <- NULL
  # } else {
  #   print("Model fitted successfully.")
  # }
  # 
  # 
  # comp_time[3] <- as.numeric(round(time_nhm,3))
  # 
  # if (!is.null(model_nhm)) {
  #   save(model_nhm, file = file.path(model_dir, "model_nhm.RData"))
  # } else {
  #   print("Model nhm is NULL; not saving.")
  # }
  # 
  # gc()
  # 
  #################
  # smms (accounts ic, accounts sm)
  #################

  # data for each observation define: patient_id, covs, time of the visit, and state
  # define the graph linking states in 'state' variable
  # specify parametric models for all transition times providing density functions as well as the corresponding survival functions
  # define starting values for opt process
  # covs number of rows equal to the number of patients and one column for each covariate
  # call: smms(startval, data, graph, covs, abs_exact = T, mc_cores = 4, hessian_matrix = T, cmethod= "" )

  # temp <- data
  # 
  # temp$state <- "dementia-free"
  # temp <- temp %>%
  #   group_by(patient_id) %>%
  #   mutate(state = ifelse(onset == 1, "dementia", state)) %>%
  #   mutate(state = ifelse(dead == 1 & row_number() == n(), "death", state)) %>%
  #   ungroup() 
  # temp <- temp[order(temp$patient_id),]
  # 
  # row_id <- temp %>%
  #   group_by(patient_id) %>%  
  #   summarise(nrows = n())
  # 
  # temp$patient_id <- rep(1:n_pats, times=as.numeric(row_id$nrows))
  # temp <- temp %>%
  #   group_by(patient_id) %>%
  #   filter(n() > 1) %>%
  #   ungroup() 
  # 
  # #kepping only needed cols and with names request by the package
  # temp <- temp %>%
  #   mutate(onset = NULL, onset_age = NULL, dead = NULL, death_time = NULL, visits = NULL)
  # temp <- temp %>%
  #   mutate(patient = patient_id, time = age, patient_id = NULL, age = NULL )
  # temp <- temp[,c(5,6,1,2,3,4)]
  # 
  # gg = graph_from_literal("dementia-free"--+"dementia"--+"death", "dementia-free"--+"death", "dementia"--+"death")
  # 
  # X_data = aggregate(temp[,c("cov1","cov2","cov3")],by=list(temp$patient),FUN=median) 
  # X_data = as.matrix(X_data[,2:4])
  # 
  # # remark: we are working with shape (a) and rate (b) parameters in which rate is in the log scale 
  # # on rate we must set covariates effect
  # # since a must be strictly positve is better to work with exp(a) and then extract shape 
  # 
  # # true model
  # f_01 = function(param, x, tt){dgompertz(tt,exp(param[1]), exp(param[2] + param[3]*x[1] + param[4]*x[2] + param[5]*x[3]))}
  # f_02 = function(param, x, tt){dgompertz(tt,exp(param[6]), exp(param[7] + param[8]*x[1] + param[9]*x[2] + param[10]*x[3]))}
  # f_12 = function(param, x, tt){dgompertz(tt,exp(param[11]), exp(param[12] + param[13]*x[1] + param[14]*x[2] + param[15]*x[3]))}
  # 
  # S_01 = function(param, x, tt){1-pgompertz(tt,exp(param[1]), exp(param[2] + param[3]*x[1] + param[4]*x[2] + param[5]*x[3]))}
  # S_02 = function(param, x, tt){1-pgompertz(tt,exp(param[6]), exp(param[7] + param[8]*x[1] + param[9]*x[2] + param[10]*x[3]))}
  # S_12 = function(param, x, tt){1-pgompertz(tt,exp(param[11]), exp(param[12] + param[13]*x[1] + param[14]*x[2] + param[15]*x[3]))}
  # 
  # # exp + covs
  # S_01 = function(param, x, t){(1-pexp(t,exp(param[1] + param[2]*x[1] + param[3]*x[2] + param[4]*x[3])))}
  # S_02 = function(param, x, t){(1-pexp(t,exp(param[5] + param[6]*x[1] + param[7]*x[2] + param[8]*x[3])))}
  # S_12 = function(param, x, t){(1-pexp(t,exp(param[9] + param[10]*x[1] + param[11]*x[2] + param[12]*x[3])))}
  # 
  # f_01 = function(param, x, t){dexp(t,exp(param[1] + param[2]*x[1] + param[3]*x[2] + param[4]*x[3]))}
  # f_02 = function(param, x, t){dexp(t, exp(param[5] + param[6]*x[1] + param[7]*x[2] + param[8]*x[3]))}
  # f_12 = function(param, x, t){dexp(t, exp(param[9] + param[10]*x[1] + param[11]*x[2] + param[12]*x[3]))}
  # 
  # # gompertz no covs
  # f_01 = function(param, x, tt){dgompertz(tt, exp(param[1]), exp(param[2]))}
  # f_02 = function(param, x, tt){dgompertz(tt, exp(param[3]), exp(param[4]))}
  # f_12 = function(param, x, tt){dgompertz(tt, exp(param[5]), exp(param[6]))}
  # 
  # S_01 = function(param, x, tt){1-pgompertz(tt, exp(param[1]), exp(param[2]))}
  # S_02 = function(param, x, tt){1-pgompertz(tt, exp(param[3]), exp(param[4]))}
  # S_12 = function(param, x, tt){1-pgompertz(tt, exp(param[5]), exp(param[6]))}
  # 
  # 
  # #just exp  
  # S_01 = function(param, x, t){(1-pexp(t,exp(param[1])))}
  # S_02 = function(param, x, t){(1-pexp(t,exp(param[2])))}
  # S_12 = function(param, x, t){(1-pexp(t,exp(param[3])))}
  # 
  # f_01 = function(param, x, t){dexp(t,exp(param[1]))}
  # f_02 = function(param, x, t){dexp(t,exp(param[2]))}
  # f_12 = function(param, x, t){dexp(t,exp(param[3]))}
  # 
  # 
  # print(names_of_survival_density(gg))
  # 
  # startval <- rep(0.5,3)
  # 
  # 
  # startval <- c(fits_gompertz[[1]]$coefficients[2], fits_gompertz[[1]]$coefficients[1], fits_gompertz[[2]]$coefficients[2], fits_gompertz[[2]]$coefficients[1], fits_gompertz[[3]]$coefficients[2], fits_gompertz[[3]]$coefficients[1])
  # model_smms <- smms_mine(startval, temp, gg, X_data, abs_exact = T, mc_cores = 4, hessian_matrix = FALSE)
  # smms(startval, temp, gg, X_data, abs_exact = T, mc_cores = 4, hessian_matrix = FALSE)

  
  ####################
  # imputation model
  ####################
  
  temp <- prepare_imputation(data, n_pats)

  
  m <- 30
  type <- "mix"
  
  error <- F
  
  time_imp <- system.time({
    tryCatch({
      results_imp<- run_imputation(temp[[1]], temp[[2]], m, type)
      avg_parameters <- results_imp[[1]]
      all_fits <- results_imp[[2]]
      
    },
    error = function(e) {
      print(paste("Error during model fitting:", e$message))
      error <<- TRUE
    })
  })[3]
  
  if (error) {
    print(paste("No imp convergence for seed:", seed))
    results_imp <- NULL
  } else {
    print("Model fitted successfully.")
  }
  
  
  comp_time[4] <- as.numeric(round(time_imp,3))
  
  if (!is.null(results_imp)) {
    save(results_imp, file = file.path(model_dir, "results_imp.RData"))
  } else {
    print("Model imp is NULL; not saving.")
  }
  
  gc()

  

  save(comp_time, file = file.path(model_dir, "computational_time.RData"))
  
  cat("models completed for seed:", seed, "\n")
  
  }
