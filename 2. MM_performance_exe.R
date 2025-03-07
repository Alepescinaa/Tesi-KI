####################################
# Upload library and data
####################################

library(fs)
library(elect)
library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(here)
library(future)
library(future.apply)
library(mstate)
library(flexsurv)
library(deSolve)
library(kableExtra)
library(hesim)

setwd(here())

load("./wrapper_MM/ground_truthMM.RData")

source_files <- c(
  "./wrapper_MM/functions_performance/check_convergence.R",
  "./wrapper_MM/functions_performance/computing_life_expectancy.R",
  "./wrapper_MM/functions_performance/compute_bias.R",
  "./wrapper_MM/functions_performance/compute_bias_rel.R",
  "./wrapper_MM/functions_performance/compute_CI.R",
  "./wrapper_MM/functions_performance/compute_coverage.R",
  "./wrapper_MM/functions_performance/compute_power.R",
  "./wrapper_MM/functions_performance/extract_comp_time.R",
  "./wrapper_MM/functions_performance/get_params_nhm.R",
  "./wrapper_MM/functions_performance/get_time.R",
  "./wrapper_MM/functions_performance/gt_flexsurv.R",
  "./wrapper_MM/functions_performance/hazards_mine.R",
  "./wrapper_MM/functions_performance/ic_comparison.R",
  "./wrapper_MM/functions_performance/ic_comparison_baseline.R",
  "./wrapper_MM/functions_performance/is.flexsurvlist.R",
  "./wrapper_MM/functions_performance/level_convergence.R",
  "./wrapper_MM/functions_performance/mean_bias_comparison.R",
  "./wrapper_MM/functions_performance/mean_bias_comparison_baseline.R",
  "./wrapper_MM/functions_performance/mean_coverage_comparison.R",
  "./wrapper_MM/functions_performance/mean_lfe_comparison.R",
  "./wrapper_MM/functions_performance/mean_power.R",
  "./wrapper_MM/functions_performance/mean_se.R",
  "./wrapper_MM/functions_performance/mean_width_ic.R",
  "./wrapper_MM/functions_performance/p.matrix.age.R",
  "./wrapper_MM/functions_performance/power_categorical.R",
  "./wrapper_MM/functions_performance/run_baseline_bias.R",
  "./wrapper_MM/functions_performance/run_performance_bias.R",
  "./wrapper_MM/functions_performance/run_performance_bias_rel.R",
  "./wrapper_MM/functions_performance/run_performance_coverage.R",
  "./wrapper_MM/functions_performance/simulation_probs.R",
  "./wrapper_MM/functions_performance/simulation_probs_nhm.R",
  "./wrapper_MM/functions_performance/standard_error.R",
  "./wrapper_MM/functions_performance/table_power.R",
  "./wrapper_MM/functions_performance/type_1_error.R",
  "./wrapper_MM/functions_performance/wrapper_convergence.R",
)

lapply(source_files, source)


# this code has to be run over each different sample size, is not taken as parameter !
# select number of patients and core to use 

n_pats <- 10000

cores <- detectCores()


###################
# load quantities #
###################

setwd(here())

if (n_pats==500){
  load("./Simulated_data_MM/simulation500_MM_all.RData")
  data <- dataset_all_MM_500
  path <- "./wrapper_MM/saved_performance_500/"
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  lapply(files, load, .GlobalEnv)
}else if (n_pats==2000){
  load("./Simulated_data_MM/simulation2K_MM_all.RData")
  data <- dataset_all_MM_2K
  path <- "./wrapper_MM/saved_performance_2K/"
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  lapply(files, load, .GlobalEnv)
}else if (n_pats==5000){
  load("./Simulated_data_MM/simulation5K_MM_all.RData")
  data <- dataset_all_MM_5K
  path <- "./wrapper_MM/saved_performance_5K/"
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  lapply(files, load, .GlobalEnv)
}else if (n_pats==10000){
  load("./Simulated_data_MM/simulation10K_MM_all.RData")
  data <- dataset_all_MM_10K
  path <- "./wrapper_MM/saved_performance_10K/"
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  lapply(files, load, .GlobalEnv)
}


# directory to save results
model_dir <- here() 
setwd(model_dir)
if (n_pats==500){
  model_dir <- paste0("wrapper_MM/saved_performance_500")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)
} else if (n_pats==2000){
  model_dir <- paste0("wrapper_MM/saved_performance_2K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)
} else if (n_pats==5000){
  model_dir <- paste0("wrapper_MM/saved_performance_5K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)  
} else if (n_pats==10000){
  model_dir <- paste0("wrapper_MM/saved_performance_10K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)  }


#######################################################
# fitting parametric model over exactly observed data #
#######################################################

# Computation of the benchmark model

plan(multisession, workers = cores)
future_lapply(1:100, function(seed) {gt_flexsurv(n_pats, seed)})

###############
# convergence #
###############

temp <- vector(mode = "list", length = 4)
convergence_schemes <- vector(mode = "list", length = 4)
hessian_schemes <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  temp[[scheme-1]] <- wrapper_convergence(n_pats, scheme, seed )
  convergence_schemes[[scheme-1]] <- temp[[scheme-1]][[1]]
  hessian_schemes[[scheme-1]] <- temp[[scheme-1]][[2]]
}

# we set combined_conv in the following way, for each method
# 0 if algorithm criteria of convergence were not met 
# 1 if the algorithm converged but not to the optimum, so no hessian exists (a bit different for coxph)
# 2 if the algorithm converged to the optimum

combined_cov <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  combined_cov[[scheme-1]] <- level_convergence(scheme)
}

setwd(here())
save(combined_cov, file = file.path(model_dir,"convergence.RData"))


######################################
# absolute bias of covariates effect #
######################################

bias_all_schemes <- vector(mode = "list", length = 4)
estimates <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  results <- data.frame(
    rate = numeric(0),
    shape = numeric(0),
    cov1 = numeric(0),
    cov2 = numeric(0),
    cov3 = numeric(0),
    `exp(cov1)` = numeric(0),
    `exp(cov2)` = numeric(0),
    `exp(cov3)` = numeric(0),
    model = character(0),
    seed = integer(0)
  )
  
  plan(multisession, workers = cores) 
  results_list <- future_lapply(1:100, function(seed) {
    temp_results <- run_performance_bias(n_pats, scheme, seed, combined_cov[[scheme - 1]])
    return(temp_results)
  })
  
  temp_bias <- list()
  temp_est <- list()
  for (seed in 1:100){
    temp_bias[[seed]] <-results_list[[seed]][[1]]
    temp_est[[seed]] <-results_list[[seed]][[2]]
  }
  
  results <- do.call(rbind, temp_bias)
  bias_all_schemes[[scheme-1]] <- results
  
  results<- do.call(rbind, temp_est)
  estimates[[scheme-1]] <- results
}

res_bias <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  res_bias[[scheme-1]] <- ic_comparison(bias_all_schemes, scheme)
}

mean_estimates <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  mean_estimates[[scheme-1]] <- mean_bias_comparison(estimates, scheme)
}

setwd(here())
save(estimates, file = file.path(model_dir,"all_estimates.RData")) 
#contains estimates of baseline parameters, covariates effect, HR for each seed

save(mean_estimates, file = file.path(model_dir,"mean_estimates.RData"))
#contains estimates of baseline parameters, covariates effect, HR averaged over the seed

save(bias_all_schemes, file = file.path(model_dir,"bias_all.RData"))
#contains absolute bias of baseline parameters, covariates effect, HR for each seed

save(res_bias, file = file.path(model_dir,"res_bias.RData"))
#contains mean and IC of absolute bias of baseline parameters and covariates effect

########################################
# relative bias of baseline parameters #
########################################

#computed over parametric models only following Gompertz distribution

baseline_all_schemes <- vector(mode = "list", length = 4)
baseline_estimates <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  results <- data.frame(
    rate = numeric(0),
    shape = numeric(0),
    model = character(0),
    seed = integer(0)
  )
  
  plan(multisession, workers = cores) 
  results_list <- future_lapply(1:100, function(seed) {
    temp_results <- run_baseline_bias(n_pats, scheme, seed, combined_cov[[scheme - 1]])
    return(temp_results)
  })
  
  temp_bias <- list()
  temp_est <- list()
  for (seed in 1:100){
    temp_bias[[seed]] <-results_list[[seed]][[1]]
    temp_est[[seed]] <-results_list[[seed]][[2]]
  }
  
  results <- do.call(rbind, temp_bias)
  baseline_all_schemes[[scheme-1]] <- results
  
  results<- do.call(rbind, temp_est)
  baseline_estimates[[scheme-1]] <- results
}

baseline_bias <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  baseline_bias[[scheme-1]] <- ic_comparison_baseline(baseline_all_schemes, scheme)
}

mean_estimates_baseline <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  mean_estimates_baseline[[scheme-1]] <- mean_bias_comparison_baseline(baseline_estimates, scheme)
}

setwd(here())
save(baseline_bias, file = file.path(model_dir,"baseline_bias.RData"))
#contains mean and IC of relative bias of baseline parameters and covariates effect

save(mean_estimates_baseline, file = file.path(model_dir,"mean_estimates_baseline.RData"))
#contains estimates of baseline parameters averaged over seeds

save(baseline_estimates, file = file.path(model_dir,"baseline_estimates_all.RData"))
#contains estimates of baseline parameters for each seed


# ######################################
#  95 % coverage of covariates effect  #
# ######################################

coverage_all_schemes <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  results_coverage <- data.frame(
    rate = numeric(0),
    shape = numeric(0),
    cov1 = numeric(0),
    cov2 = numeric(0),
    cov3 = numeric(0),
    model = character(0),
    seed = integer(0)
  )

  plan(multisession, workers = cores)
  results_list <- future_lapply(1:100, function(seed) {
    temp_results <- run_performance_coverage(n_pats, scheme, seed, combined_cov[[scheme-1]])
    return(temp_results)
  })

  results_coverage <- do.call(rbind, results_list)

  coverage_all_schemes[[scheme-1]] <- results_coverage

}

res_cov <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  res_cov[[scheme-1]] <- mean_coverage_comparison(coverage_all_schemes, scheme)
}

res_cov_ic <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  res_cov_ic[[scheme-1]] <- ic_comparison(coverage_all_schemes, scheme)
}

setwd(here())
save(coverage_all_schemes, file = file.path(model_dir,"all_coverage.RData"))
#contains coverage of baseline parameters and covariates effects for each seed

save(res_cov_ic, file = file.path(model_dir,"95%coverage.RData"))
#contains coverage mean and IC of baseline parameters and covariates effects 


########################################
# standard error of covariates effect  #
########################################

se_all_schemes <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  plan(multisession, workers = cores) 
  results_list <- future_lapply(1:100, function(seed) {
    temp_results <- standard_error(n_pats, scheme, seed, combined_cov[[scheme-1]])
    return(temp_results)
  })
  
  results_se<- do.call(rbind, results_list)
  
  se_all_schemes[[scheme-1]] <- results_se
  
}

se_mean<- vector(mode = "list", length = 4)
for (scheme in 2:5){
  se_mean[[scheme-1]] <- mean_standard_error(se_all_schemes, scheme)
}

setwd(here())
save(se_all_schemes, file = file.path(model_dir,"se_all.RData"))
#contains standard error of covariates effects for each seed

save(se_mean, file = file.path(model_dir,"mean_se.RData"))
#contains mean standard error of covariates effects 


####################
# life expectancy #
####################

lfe_bias <- vector(mode = "list", length = 4)
lfe_estimates <- vector(mode = "list", length = 4)
covs <-  data.frame(
  cov1 = 0.35, 
  cov2 = 0.15,
  cov3 = 0.30
)

for (scheme in 2:5){
  results <- data.frame(
    lfe = numeric(0),
    model = character(0),
    seed = integer(0)
  )
  # I've been using msm age without accounting for age as time varying
  plan(multisession, workers = cores) 
  results_list <- future_lapply(1:100, function(seed) {
    t_start <- 60
    temp_results <- computing_life_expectancy(n_pats, scheme, seed, combined_cov[[scheme-1]], t_start, covs)
    return(temp_results)
  })
  
  temp_bias <- list()
  temp_est <- list()
  for (seed in 1:100){
    temp_est[[seed]] <-results_list[[seed]][[1]]
    temp_bias[[seed]] <-results_list[[seed]][[2]]
  }
  
  results <- do.call(rbind, temp_bias)
  lfe_bias[[scheme-1]] <- results
  
  results<- do.call(rbind, temp_est)
  lfe_estimates[[scheme-1]] <- results
}

res_bias_lfe <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  res_bias_lfe[[scheme-1]] <- mean_lfe_comparison(lfe_bias, scheme)
}

mean_estimates_lfe<- vector(mode = "list", length = 4)
for (scheme in 2:5){
  mean_estimates_lfe[[scheme-1]] <- mean_lfe_comparison(lfe_estimates, scheme)
}

setwd(here())
#save(gt_tls, file = file.path(model_dir,"gt_tls.RData"))
# contains the ground truth averaged time spent in each state 

save(lfe_estimates, file = file.path(model_dir,"all_estimates_lfe.RData")) 
# contains average time spent in each state for each seed
save(mean_estimates_lfe, file = file.path(model_dir,"mean_estimates_lfe.RData"))
# contains average time spent in each state averaged over seeds

save(res_bias_lfe, file = file.path(model_dir,"bias_lfe.RData"))
# contains average bias of the time spent in each state averaged over seeds

##########################
# Power and type 1 error #
##########################

# Type 1 error -> refuse H0|H0 true -> significant|no significant
# Power -> refuse H0|H0 false -> significant|significant

significancy <- vector(mode = "list", length = 4)
significancy_all <- vector(mode = "list", length = 4)

alpha <- 0.05

for (scheme in 2:5){
  plan(multisession, workers = cores) 
  results_list <- future_lapply(1:100, function(seed) {
    temp_results <- compute_power(n_pats, scheme, seed, combined_cov[[scheme - 1]], alpha)
    return(temp_results)
  })
  
  results<- do.call(rbind, results_list)
  significancy_all[[scheme-1]] <- results
  significancy[[scheme-1]] <- mean_power(results, scheme)
  

}

significant_covs <- data.frame("cov1"= c(0,1,1), "cov2"= c(1,1,0), "cov3"=c(1,1,1), "transition"=c(1,2,3))

setwd(here())
save(significancy, file = file.path(model_dir,"significancy.RData"))
# mean probability that the covariates are predicted to be significant 

save(significancy_all, file = file.path(model_dir,"significancy_all.RData"))
# whether each covariate is predicted to be significant for each seed

save(significant_covs, file = file.path(model_dir,"significancy_covs.RData"))
# presence/absence of covariates effect in the simulation process


pw2 <- power_categorical(significant_covs, significancy, scheme=2)
pw3 <- power_categorical(significant_covs, significancy, scheme=3)
pw4 <- power_categorical(significant_covs, significancy, scheme=4)
pw5 <- power_categorical(significant_covs, significancy, scheme=5)

setwd(model_dir)
plots_power <- list(pw2,pw3,pw4,pw5)
save(plots_power, file = "power.RData" )
# contains mean power for each transition


setwd(here())

err2 <- type_1_error(significant_covs, significancy, scheme=2)
err3 <- type_1_error(significant_covs, significancy, scheme=3)
err4 <- type_1_error(significant_covs, significancy, scheme=4)
err5 <- type_1_error(significant_covs, significancy, scheme=5)

setwd(model_dir)
plots_errorI <- list(err2,err3,err4,err5)
save(plots_errorI, file = "typeIerr.RData" )
#contains mean type I error for each feasible transition

setwd(here())

######################
# computational time #
######################

ct_all_schemes <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  plan(multisession, workers = cores) 
  results_list <- future_lapply(1:100, function(seed) {
    temp_results <- extract_comp_time(n_pats, scheme, seed)
    return(temp_results)
  })
  
  results_ct <- do.call(rbind, results_list)
  
  ct_all_schemes[[scheme-1]] <- results_ct
}

setwd(here())
save(ct_all_schemes, file = file.path(model_dir,"comp_time.RData"))
#contains mean computational time 


##########
# Plots
##########

dir <- here()
dir <- here("wrapper_MM", "Plots")
setwd(dir)

titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
plot1 <- plot_convergence(2, titles)
plot2 <- plot_convergence(3, titles)
plot3 <- plot_convergence(4, titles)
plot4 <- plot_convergence(5, titles)

combined_plot <- (plot1 + theme(legend.position = "none")) +
  (plot2 + theme(legend.position = "none")) +
  (plot3 + theme(legend.position = "none")) +
  (plot4 + theme(legend.position = "none")) +              
  plot_layout(guides = "collect")  

combined_plot <- combined_plot & theme(legend.position = "right")

ggsave("convM.png", plot = combined_plot, path = NULL, width = 10, height = 7) 

 
titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
pb2 <- plot_bias(2, titles)
plot_bias(3, titles)
pb4 <- plot_bias(4, titles)
plot_bias(5, titles)

ggsave("bias2.png", plot = pb2, path = NULL, width = 9, height = 7)
ggsave("bias4.png", plot = pb4, path = NULL, width = 9, height = 7)

# for (scheme in 2:5) {
#   for (transition in 1:3) {
#     print(plot_bias_rel(scheme, titles, transition))
#   }
# }

# scheme <- 4 #2:5
# list_data <- prepare_data_boxplot(scheme)
# 
# parameters <- c("cov1", "cov2", "cov3", "rate", "shape")
# plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[1])
# plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[2])
# plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[3])
# plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[4])
# plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[5])

titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
cov2 <- plot_coverage(2, titles)
plot_coverage(3, titles)
cov4 <- plot_coverage(4, titles)
plot_coverage(5, titles)

ggsave("cov2.png", plot = cov2, path = NULL, width = 9, height = 7)
ggsave("cov4.png", plot = cov4, path = NULL, width = 9, height = 7)

titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
plot_se(se_mean,2, titles)
plot_se(se_mean,3, titles)
plot_se(se_mean,4, titles)
plot_se(se_mean,5, titles)

ggsave("width4.png", plot = w4, path = NULL, width = 9, height = 7)

plfe <- plot_lfe(0)
plfe_dem <- plot_lfe(1) 

#relative bias
plfe_b <- plot_lfe_bias(0)
plfe_dem_b <- plot_lfe_bias(1)

dir <- here()
dir <- here("wrapper_MM", "Plots")
setwd(dir)

ggsave("lfeM.png", plot = plfe, path = NULL, width = 10, height = 7) 
ggsave("lfe_demM.png", plot = plfe_dem, path = NULL, width = 10, height = 7) 
ggsave("lfeM_b.png", plot = plfe_b, path = NULL, width = 10, height = 7) 
ggsave("lfe_demM_b.png", plot = plfe_dem_b, path = NULL, width = 10, height = 7) 



# green type one error
# yellow power
table_power(significant_covs, significancy, scheme=2)
table_power(significant_covs, significancy, scheme=3)
t4 <- table_power(significant_covs, significancy, scheme=4)
table_power(significant_covs, significancy, scheme=5)

save_kable(t4, file = "power.html")
titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
plot1 <- plot_ct(2, titles, combined_cov[[1]])
plot2 <- plot_ct(3, titles, combined_cov[[2]])
plot3 <- plot_ct(4, titles, combined_cov[[3]])
plot4 <- plot_ct(5, titles, combined_cov[[4]])


ct_plot <- (plot1 + theme(legend.position = "none")) +
  (plot2 + theme(legend.position = "none")) +
  (plot3 + theme(legend.position = "none")) +
  (plot4 + theme(legend.position = "none")) +       
  plot_layout(guides = "collect")  

ct_plot <- ct_plot & theme(legend.position = "left")

ggsave("ct.png", plot = ct_plot, path = NULL, width = 9, height = 7)


