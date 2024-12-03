####################################
# Upload library and data
####################################

library(fs)
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

setwd(here())

load("./wrapper_SM/ground_truthSM.RData")
# ground_truth_params <- ground_truth_params[,c(2,1,3,4,5)]
# ground_truth_params <- cbind(ground_truth_params, exp(ground_truth_params[,3]), exp(ground_truth_params[,4]),exp(ground_truth_params[,5]))
# colnames(ground_truth_params) <- c("rate","shape","cov1","cov2","cov3", "exp(cov1)", "exp(cov2)", "exp(cov3)")
# rownames(ground_truth_params) <- c("1->2", "1->3", "2->3")
# model_dir <-"/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi-KI/wrapper_SM"
# save(ground_truth_params, file = file.path(model_dir,"ground_truthSM.RData"))

source_files <- c(
  "./wrapper_SM/functions_performance/gt_flexsurv.R",
  "./wrapper_SM/functions_performance/compute_bias.R",
  "./wrapper_SM/functions_performance/run_performance_bias.R",
  "./wrapper_SM/functions_performance/mean_bias_comparison.R",
  "./wrapper_SM/functions_performance/wrapper_convergence.R",
  "./wrapper_SM/functions_performance/check_convergence.R",
  "./wrapper_SM/functions_performance/level_convergence.R",
  "./wrapper_SM/functions_performance/compute_coverage.R",
  "./wrapper_SM/functions_performance/compute_CI.R",
  "./wrapper_SM/functions_performance/run_performance_coverage.R",
  "./wrapper_SM/functions_performance/get_params_nhm.R",
  "./wrapper_SM/functions_performance/mean_coverage_comparison.R",
  "./wrapper_SM/functions_performance/computing_life_expectancy.R",
  "./wrapper_SM/functions_performance/mean_lfe_comparison.R",
  "./wrapper_SM/functions_performance/extract_comp_time.R",
  "./wrapper_MM/functions_performance/width_ic.R",
  "./wrapper_MM/functions_performance/mean_width_ic.R",
  "./wrapper_MM/functions_performance/compute_power.R",
  "./wrapper_MM/functions_performance/mean_power.R",
  "./wrapper_MM/functions_performance/table_power.R",
 "./wrapper_SM/functions_performance/plot_convergence.R",
 "./wrapper_SM/functions_performance/plot_bias.R",
"./wrapper_SM/functions_performance/ic_comparison.R",
"./wrapper_SM/functions_performance/plot_coverage.R",
"./wrapper_SM/functions_performance/plot_width.R",
"./wrapper_SM/functions_performance/plot_ct.R"
)
  
  # "./wrapper_SM/functions_performance/compute_bias_rel.R",
  # "./wrapper_SM/functions_performance/run_performance_bias_rel.R",
  
  # "./wrapper_SM/functions_performance/plot_bias.R",
  # "./wrapper_SM/functions_performance/plot_bias_rel.R",
  # 
  # "./wrapper_SM/functions_performance/plot_boxplot.R",
  # "./wrapper_SM/functions_performance/plot_ct.R",
  # "./wrapper_SM/functions_performance/prepare_data_boxplot.R",
  # "./wrapper_SM/functions_performance/ic_comparison.R",
  # "./wrapper_SM/functions_performance/is.flexsurvlist.R",
  # "./wrapper_SM/functions_performance/plot_bias_lfe.R",
  # 



lapply(source_files, source)


# this code has to be run over each different sample size, is not taken as parameter !
# select number of patients and core to use 

n_pats <- 5000
cores <- 4


###################
# load quantities #
###################

setwd(here())

if (n_pats==500){
  load("./Simulated_data_SM/simulation500_SM_all.RData")
  data <- dataset_all_SM_500
  path <- "./wrapper_SM/saved_performance_500/"
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  lapply(files, load, .GlobalEnv)
}else if (n_pats==2000){
  load("./Simulated_data_SM/simulation2K_SM_all.RData")
  data <- dataset_all_SM_2K
  path <- "./wrapper_SM/saved_performance_2K/"
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  lapply(files, load, .GlobalEnv)
}else if (n_pats==5000){
  load("./Simulated_data_SM/simulation5K_SM_all.RData")
  data <- dataset_all_SM_5K
  path <- "./wrapper_SM/saved_performance_5K/"
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  lapply(files, load, .GlobalEnv)
}else if (n_pats==10000){
  load("./Simulated_data_SM/simulation10K_SM_all.RData")
  data <- dataset_all_SM_10K
  path <- "./wrapper_SM/saved_performance_10K/"
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  lapply(files, load, .GlobalEnv)
}

# directory to save things
model_dir <- here()
setwd(model_dir)
if (n_pats==500){
  model_dir <- paste0("wrapper_SM/saved_performance_500")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)
} else if (n_pats==2000){
  model_dir <- paste0("wrapper_SM/saved_performance_2K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)
} else if (n_pats==5000){
  model_dir <- paste0("wrapper_SM/saved_performance_5K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)  
} else if (n_pats==10000){
  model_dir <- paste0("wrapper_SM/saved_performance_10K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)  }

#######################################################
# fitting parametric model over exactly observed data #
#######################################################

# plan(multisession, workers = cores)
# future_lapply(1:100, function(seed) {gt_flexsurv(n_pats, seed)})
# for (seed in 1:100){
#   gt_flexsurv(n_pats, seed)
#   print(seed)
# }

########################
# check of convergence #
########################

temp <- vector(mode = "list", length = 4)
convergence_schemes <- vector(mode = "list", length = 4)
hessian_schemes <- vector(mode = "list", length = 4)

for (scheme in 2:3){
  temp[[scheme-1]] <- wrapper_convergence(n_pats, scheme, seed ) 
  convergence_schemes[[scheme-1]] <- temp[[scheme-1]][[1]]
  hessian_schemes[[scheme-1]] <- temp[[scheme-1]][[2]]
}

# we set combined_conv in the following way, for each method
# 0 if algorithm criteria of convergence were not met 
# 1 if the algorithm converged but not to the optimum, so no hessian exists (a bit different for coxph)
# 2 if the algorithm converged to the optimum

combined_cov <- vector(mode = "list", length = 4)

for (scheme in 2:3){
  combined_cov[[scheme-1]] <- level_convergence(scheme)
}

 
save(combined_cov, file = file.path(model_dir,"convergence.RData"))

###################
# bias comparison #
###################

bias_all_schemes <- vector(mode = "list", length = 4)
estimates <- vector(mode = "list", length = 4)

for (scheme in 2:3){
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
for (scheme in 2:3){
  res_bias[[scheme-1]] <- mean_bias_comparison(bias_all_schemes, scheme)
}

mean_estimates <- vector(mode = "list", length = 4)
for (scheme in 2:3){
  mean_estimates[[scheme-1]] <- mean_bias_comparison(estimates, scheme)
}

save(mean_estimates, file = file.path(model_dir,"mean_estimates.RData"))
save(bias_all_schemes, file = file.path(model_dir,"bias_all.RData"))
save(res_bias, file = file.path(model_dir,"res_bias.RData"))


############################
# relative bias comparison #
############################

# rel_bias_all_schemes <- vector(mode = "list", length = 4)
# 
# for (scheme in 2:5){
#   results_bias_rel <- data.frame(
#     rate = numeric(0),
#     shape = numeric(0),
#     cov1 = numeric(0),
#     cov2 = numeric(0),
#     cov3 = numeric(0),
#     `exp(cov1)` = numeric(0),
#     `exp(cov2)` = numeric(0),
#     `exp(cov3)` = numeric(0),
#     model = character(0),
#     seed = integer(0)
#   )
#   
#   plan(multisession, workers = cores) 
#   results_list <- future_lapply(1:100, function(seed) {
#     temp_results <- run_performance_bias_rel(n_pats, scheme, seed, combined_cov[[scheme - 1]])
#     return(temp_results)
#   })
#   
#   results_bias_rel <- do.call(rbind, results_list)
#   
#   rel_bias_all_schemes[[scheme-1]] <- results_bias_rel
# }
# 
# res_bias_rel <- vector(mode = "list", length = 4)
# for (scheme in 2:5){
#   res_bias_rel[[scheme-1]] <- mean_bias_comparison(rel_bias_all_schemes, scheme)
# }
# 
# save(rel_bias_all_schemes, file = file.path(model_dir,"bias_all_rel.RData"))
# save(res_bias_rel, file = file.path(model_dir,"res_bias_rel.RData"))

###########################
# 95% coverage comparison #
###########################

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

save(coverage_all_schemes, file = file.path(model_dir,"all_coverage.RData"))
save(res_cov_ic, file = file.path(model_dir,"95%coverage.RData"))

###############################
#  width confidence intervals #
###############################

width_all_schemes <- vector(mode = "list", length = 4)

for (scheme in 2:5){
  plan(multisession, workers = cores) 
  results_list <- future_lapply(1:100, function(seed) {
    temp_results <- width_ic(n_pats, scheme, seed, combined_cov[[scheme-1]])
    return(temp_results)
  })
  
  results_width<- do.call(rbind, results_list)
  
  width_all_schemes[[scheme-1]] <- results_width
  
}

mean_width <- vector(mode = "list", length = 4)
for (scheme in 2:5){
  mean_width[[scheme-1]] <- mean_width_ic(width_all_schemes, scheme)
}

save(width_all_schemes, file = file.path(model_dir,"width_all_ic.RData"))
save(mean_width, file = file.path(model_dir,"width_ic.RData"))

###################
# life expectancy #
###################

cores_lfe=4
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
  
  plan(multisession, workers = cores) 
  results_list <- future_lapply(1:10, function(seed) {
    temp <- data[[seed]][[scheme]]
    t_start <- min(temp$age)
    temp_results <- computing_life_expectancy(n_pats, scheme, seed, combined_cov[[scheme-1]], t_start, covs)
    return(temp_results)
  })
  
  temp_bias <- list()
  temp_est <- list()
  for (seed in 1:10){
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

save(mean_estimates_lfe, file = file.path(model_dir,"mean_estimates_lfe.RData"))
save(res_bias_lfe, file = file.path(model_dir,"bias_lfe.RData"))

##########################
# Power and type 1 error #
##########################

# In the table I will represent P(p_i<=alpha) i.e. percentage of times for which the covs effect is significant
# the yellow values represent when the covariate was significant in the ground truth 
# H0 is variable not significant
# Type 1 error -> refuse H0|H0 true -> significant|no sign
# Type 2 error ->  accept H0| H0 false -> no sign|sign
# Power -> refuse H0|H0 false -> significant|sign

significancy_all <- vector(mode = "list", length = 4)
significancy <- vector(mode = "list", length = 4)
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

save(significancy, file = file.path(model_dir,"significancy.RData"))
save(significancy_all, file = file.path(model_dir,"significancy_all.RData"))

significant_covs <- data.frame("cov1"= c(0,1,1), "cov2"= c(1,1,0), "cov3"=c(1,1,1), "transition"=c(1,2,3))

# green type one error
# yellow power
table_power(significant_covs, significancy, scheme=2)
table_power(significant_covs, significancy, scheme=3)
table_power(significant_covs, significancy, scheme=4)
table_power(significant_covs, significancy, scheme=5)



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

save(ct_all_schemes, file = file.path(model_dir,"comp_time.RData"))


#########
# Plots #
#########

titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
plot1 <- plot_convergence(2, titles)
plot2 <- plot_convergence(3, titles)
plot3 <- plot_convergence(4, titles)
plot4 <- plot_convergence(5, titles)

plot1 + plot2 + plot3 + plot4 + 
  plot_layout(ncol = 2, nrow = 2, 
              widths = c(0.6, 0.6), 
              heights = c(0.8, 0.8)) +
  theme(plot.margin = margin(10, 10, 10, 10))

titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
plot_bias(2, titles)
plot_bias(3, titles)
plot_bias(4, titles)
plot_bias(5, titles)

for (scheme in 2:5) {
  for (transition in 1:3) {
    print(plot_bias_rel(scheme, titles, transition))
  }
}

scheme <- 4 #2:5
list_data <- prepare_data_boxplot(scheme)

parameters <- c("cov1", "cov2", "cov3", "rate", "shape")
plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[1])
plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[2])
plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[3])
plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[4])
plot_boxplot(list_data[[1]], list_data[[2]], list_data[[3]], parameters[5])

titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
plot_coverage(2, titles)
plot_coverage(3, titles)
plot_coverage(4, titles)
plot_coverage(5, titles)

titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
plot_width(mean_width,2, titles)
plot_width(mean_width,3, titles)
plot_width(mean_width,4, titles)
plot_width(mean_width,5, titles)

titles <-c("(1 year)", "(3 years)", "(3-6 years)", "EHR")
plot_bias_lfe(res_bias_lfe, 2, titles) 
plot_bias_lfe(res_bias_lfe, 3, titles)
plot_bias_lfe(res_bias_lfe, 4, titles)
plot_bias_lfe(res_bias_lfe, 5, titles)

titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
plot1 <- plot_ct(2, titles, combined_cov[[1]])
plot2 <- plot_ct(3, titles, combined_cov[[2]])
plot3 <- plot_ct(4, titles, combined_cov[[3]])
plot4 <- plot_ct(5, titles, combined_cov[[4]])
plot1+plot2+plot3+plot4

dir <- here()
dir <- paste0("wrapper_SM/Plots_500")
ggsave("plot3.png", path=dir, width=5, height=8, bg = "transparent")

# keep in mind that these estimates of the bias are accounted only for the models for which convergence at optimum is met
