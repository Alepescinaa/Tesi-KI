wrapper_functions_MM <- function(data,n_pats,seed,cores_nhm){
  
  comp_time <- numeric()
  setwd(here())

  if (n_pats==500){
    if (scheme==2){
      model_dir <- paste0("wrapper_MM/results_500/saved_models_scheme2/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==3){
      model_dir <- paste0("wrapper_MM/results_500/saved_models_scheme3/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==4){
      model_dir <- paste0("wrapper_MM/results_500/saved_models_scheme4/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  
    } else if (scheme==5){
      model_dir <- paste0("wrapper_MM/results_500/saved_models_scheme5/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  }
  } else if (n_pats==2000){
    if (scheme==2){
      model_dir <- paste0("wrapper_MM/results_2K/saved_models_scheme2/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==3){
      model_dir <- paste0("wrapper_MM/results_2K/saved_models_scheme3/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==4){
      model_dir <- paste0("wrapper_MM/results_2K/saved_models_scheme4/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  
    } else if (scheme==5){
      model_dir <- paste0("wrapper_MM/results_2K/saved_models_scheme5/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  }
  } else if (n_pats==5000){
    if (scheme==2){
      model_dir <- paste0("wrapper_MM/results_5K/saved_models_scheme2/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==3){
      model_dir <- paste0("wrapper_MM/results_5K/saved_models_scheme3/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==4){
      model_dir <- paste0("wrapper_MM/results_5K/saved_models_scheme4/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  
    } else if (scheme==5){
      model_dir <- paste0("wrapper_MM/results_5K/saved_models_scheme5/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  }
  } else if (n_pats==10000){
    if (scheme==2){
      model_dir <- paste0("wrapper_MM/results_10K/saved_models_scheme2/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==3){
      model_dir <- paste0("wrapper_MM/results_10K/saved_models_scheme3/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)
    } else if (scheme==4){
      model_dir <- paste0("wrapper_MM/results_10K/saved_models_scheme4/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  
    } else if (scheme==5){
      model_dir <- paste0("wrapper_MM/results_10K/saved_models_scheme5/seed_", seed)
      dir.create(model_dir, showWarnings = FALSE, recursive= T)  }
  }

  n_pats <- length(unique(data$patient_id))
  
  #####################
  # coxph model
  #####################
  
  temp <-  prepare_coxph_flex(data, n_pats)
  temp <- expand.covs(temp,c("cov1","cov2", "cov3"))


  time_cox<- system.time({
  model_cox <- coxph(Surv(Tstart,Tstop,status) ~ cov1.1 + cov2.1 + cov3.1 + cov1.2 + cov2.2 + cov3.2 + cov1.3 + cov2.3 + cov3.3 + strata(trans), data = temp, method="breslow")
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

  gc()

  ######################
  # msm model
  ######################

  temp <- prepare_msm(data)

  Q <- rbind(c(0, 1, 1),
             c(0, 0, 1),
             c(0, 0, 0))

  time_msm<- system.time({
    model.msm <- msm(state ~ age,
                     subject = patient_id,
                     data = temp,
                     qmatrix = Q,
                     covariates = ~ cov1 + cov2 + cov3 ,
                     gen.inits = TRUE,
                     censor = 99,
                     control = list(fnscale = 1000, maxit = 1000),
                     deathexact = TRUE)
  })[3]

  comp_time[3] <- as.numeric(round(time_msm,3))

  initial_guess_age <- qmatrix.msm(model.msm)$estimates
  msm_estimates <- model.msm$estimates.t

  save(model.msm, file = file.path(model_dir, "msm_model.RData"))

  #####################
  # msm + age model
  #####################

  temp <- prepare_msm(data)



  time_msm_age<- system.time({
    model.msm_age <- msm(state ~ age,
                         subject = patient_id,
                         data = temp,
                         qmatrix = initial_guess_age,
                         covariates = ~ cov1 + cov2 + cov3 + age ,
                         gen.inits= FALSE,
                         control = list(fnscale = 1000, maxit = 1000),
                         censor = 99,
                         center= FALSE,
                         deathexact = TRUE)
  })[3]

  comp_time[4] <- as.numeric(round(time_msm_age,3))

  save(model.msm_age, file = file.path(model_dir, "model_msm_age.RData"))

  gc()

  ######################
  # nhm model
  ######################

  temp <- prepare_msm(data)

  tmat_1 <- rbind(c(0,1,2),c(0,0,3),rep(0,3))
  tmat_2 <- rbind(c(0,4,5),c(0,0,6),rep(0,3))
  tmat_3 <- rbind(c(0,7,8),c(0,0,9),rep(0,3))

  temp$patient_id <- as.factor(temp$patient_id)
  temp=as.data.frame(temp)

  find_splits <- function(age) {
    quantiles <- quantile(age, probs = seq(0, 1, 0.05))
    return(quantiles[-c(1, length(quantiles))])
  }

  split_points <- find_splits(temp$age)[10:19]

  initial_guess <-  append(msm_estimates, rep(0.5, 3), after = 3)
  # we have estimates for rate and covs, so we add initial estimate to 0.5 of shape -> not helping
  error <- F

  time_nhm <- system.time({
    tryCatch({
      object_nhm <- model.nhm(state ~ age,
                              subject = patient_id,
                              data = temp,
                              trans = tmat_1,
                              nonh = tmat_1,
                              type = "gompertz",
                              covariates = c("cov1", "cov2", "cov3"),
                              covm = list(cov1= tmat_1, cov2=tmat_2, cov3=tmat_3),
                              censor.states = c(1,2),
                              censor = 99,
                              death = T,
                              death.states = c(3))


      model_nhm <- nhm(object_nhm,
                       gen_inits = TRUE,
                       #initial = initial_guess,
                       score_test = FALSE,
                       control = nhm.control(ncores = cores_nhm, obsinfo = FALSE, coarsen = T, coarsen.vars = c(1), coarsen.lv = 10, splits = split_points))
    },
    error = function(e) {
      print(paste("Error during model fitting:", e$message))
      error <<- TRUE
    })
  })[3]

  if (error) {
    print(paste("No convergence for seed:", seed))
    model_nhm <- NULL
  } else {
    print("Model fitted successfully.")
  }


  comp_time[5] <- as.numeric(round(time_nhm,3))

  if (!is.null(model_nhm)) {
    save(model_nhm, file = file.path(model_dir, "model_nhm.RData"))
  } else {
    print("Model is NULL; not saving.")
  }

  gc()


  ####################
  # imputation model
  ####################

  temp <- prepare_imputation(data, n_pats)
  m <- 30
  type <- "forward"

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
    print(paste("No convergence for seed:", seed))
    results_imp <- NULL
  } else {
    print("Model fitted successfully.")
  }

  comp_time[6] <- as.numeric(round(time_imp,3))

  save(results_imp, file = file.path(model_dir, "results_imp.RData"))

  gc()

  ########################
  # Smoothhazard
  #########################

  # temp <- prepare_SmoootHazard(data, n_pats)
  #
  # error <- FALSE
  #
  # time_smh <- system.time({
  #   tryCatch({
  #     model_smootHazard <- idm(
  #       formula01 = Hist(time = list(l, r), event = onset, entry = age) ~ cov1+cov2+cov3,
  #       formula02 = Hist(time = death_time, event = dead, entry = age) ~ cov1+cov2+cov3,
  #       method = "Splines",
  #       data = temp,
  #       CV = TRUE,
  #       maxiter = 500,
  #       print.iter = TRUE
  #     )
  #   }, error = function(e) {
  #     error <<- TRUE  # Set the error flag
  #     message("An error occurred: ", e$message)
  #   })
  # })[3]
  #
  #
  # if (error) {
  #   print("The model fitting encountered an error.")
  # }
  #
  # comp_time[7] <- as.numeric(round(time_smh,3))
  #
  # if (!is.null(model_smootHazard)) {
  #   save(model_smootHazard, file = file.path(model_dir, "model_smootHazard.RData"))
  # } else {
  #   print("Model is NULL; not saving.")
  # }
  #
  #
   save(comp_time, file = file.path(model_dir, "computational_time.RData"))
  #
  cat("models completed for seed:", seed, "\n")
  
  }
