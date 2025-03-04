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
library(kableExtra)
library(ggridges)
library(scales)

setwd(here())

load("./wrapper_MM/ground_truthMM.RData")

source_files <- c(
  "./wrapper_MM/functions_performance/plot_convergence.R",
  "./wrapper_MM/functions_performance/plot_bias.R",
  "./wrapper_MM/functions_performance/plot_bias_rel.R",
  "./wrapper_MM/functions_performance/plot_coverage.R",
  "./wrapper_MM/functions_performance/plot_boxplot.R",
  "./wrapper_MM/functions_performance/plot_ct.R",
  "./wrapper_MM/functions_performance/table_power.R",
  "./wrapper_MM/functions_performance/width_ic.R",
  "./wrapper_MM/functions_performance/plot_width.R",
  "./wrapper_MM/functions_performance/plot_lfe.R",
  "./wrapper_MM/functions_performance/plot_lfe_bias.R",
  "./wrapper_MM/functions_performance/power_categorical.R",
  "./wrapper_MM/functions_performance/type_1_error.R",
  "./wrapper_MM/functions_performance/plot_distribution.R"
  
)


lapply(source_files, source)


n_pats <- 500
cores <- 4

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

###########################
# directory to save plots #
###########################

model_dir <- here() 
setwd(model_dir)
if (n_pats==500){
  model_dir <- paste0("wrapper_MM/plots_500")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)
} else if (n_pats==2000){
  model_dir <- paste0("wrapper_MM/plots_2K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)
} else if (n_pats==5000){
  model_dir <- paste0("wrapper_MM/plots_5K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)  
} else if (n_pats==10000){
  model_dir <- paste0("wrapper_MM/plots_10K")
  dir.create(model_dir, showWarnings = FALSE, recursive= T)  }




##########
# Plots
##########

tot_convergence <- rbind(combined_cov[[1]],combined_cov[[2]],combined_cov[[3]],combined_cov[[4]])
tot_convergence$scheme <- rep(2:5,each=100)
tot_convergence$seed <- rep(1:100,times=4)
tot_convergence <- tot_convergence %>%
  mutate(a=coxph, b=flexsurv, c=msm, d=msm_age, e=nhm, f=imputation, coxph=NULL, flexsurv=NULL, msm=NULL, msm_age=NULL, nhm=NULL, imputation=NULL)
conv_long <- tot_convergence %>%
  pivot_longer(cols = c("a","b","c","d","e","f"),
               names_to = "model",
               values_to = "convergence")
conv_summary <- conv_long %>%
  group_by(scheme, model, convergence) %>%
  summarise(count = n(), .groups = "drop")%>%
  group_by(scheme, model) %>%
  mutate(percentage = count / sum(count) * 100)

conv_table <- conv_summary %>%
  select(scheme, model, convergence, percentage) %>%
  pivot_wider(names_from = convergence, values_from = percentage, values_fill = 0) %>%
  rename( `Convergence 1` = `1`, `Convergence 2` = `2`)

setwd(model_dir)
save(conv_table, file = "convergence_10K.RData" )

# titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
# plot1 <- plot_convergence(2, titles)
# plot2 <- plot_convergence(3, titles)
# plot3 <- plot_convergence(4, titles)
# plot4 <- plot_convergence(5, titles)
# 
# combined_plot <- (plot1 + theme(legend.position = "none")) +
#   (plot2 + theme(legend.position = "none")) +
#   (plot3 + theme(legend.position = "none")) +
#   (plot4 + theme(legend.position = "none")) +              
#   plot_layout(guides = "collect")  
# 
# combined_plot <- combined_plot & theme(legend.position = "right")
# 
# ggsave("convergence.png", plot = combined_plot, path = model_dir, width = 10, height = 7) 


titles <- c("1 year observation scheme", "3 years observation scheme", "3-6 years observation scheme", "irregular observation scheme")
pb2 <- plot_bias(res_bias, 2, titles)
pb3 <- plot_bias(res_bias, 3, titles)
pb4 <- plot_bias(res_bias, 4, titles)
pb5 <- plot_bias(res_bias, 5, titles)

titles <- c("1 year observation scheme", "3 years observation scheme", "3-6 years observation scheme", "irregular observation scheme")
pb2 <- plot_bias_baseline(baseline_bias, 2, titles)
pb3 <- plot_bias_baseline(baseline_bias, 3, titles)
pb4 <- plot_bias_baseline(baseline_bias, 4, titles)
pb5 <- plot_bias_baseline(baseline_bias, 5, titles)

setwd(here())
ggsave("bias2_cov.png", plot = pb2[[1]], path = model_dir, width = 9, height = 7)
ggsave("bias2_base.png", plot = pb2[[2]], path = model_dir, width = 9, height = 7)
ggsave("bias3_cov.png", plot = pb3[[1]], path = model_dir, width = 9, height = 7)
ggsave("bias3_base.png", plot = pb3[[2]], path = model_dir, width = 9, height = 7)
ggsave("bias4_cov.png", plot = pb4[[1]], path = model_dir, width = 9, height = 7)
ggsave("bias4_base.png", plot = pb4[[2]], path = model_dir, width = 9, height = 7)
ggsave("bias5_cov.png", plot = pb5[[1]], path = model_dir, width = 9, height = 7)
ggsave("bias5_base.png", plot = pb5[[2]], path = model_dir, width = 9, height = 7)

plots_bias <- list( pb2[[1]], pb2[[2]], pb3[[1]], pb3[[2]], pb4[[1]], pb4[[2]], pb5[[1]], pb5[[2]])
setwd(model_dir)
save(plots_bias, file = "bias_500.RData" )

plots_bias <- list( pb2,  pb3, pb4, pb5)
setwd(model_dir)
save(plots_bias, file = "bias_baseline_5000.RData" )


cov2 <- plot_coverage(2, titles)
cov3 <- plot_coverage(3, titles)
cov4 <- plot_coverage(4, titles)
cov5 <- plot_coverage(5, titles)

setwd(here())
ggsave("cov2.png", plot = cov2, path = model_dir, width = 9, height = 7)
ggsave("cov3.png", plot = cov3, path = model_dir, width = 9, height = 7)
ggsave("cov4.png", plot = cov4, path = model_dir, width = 9, height = 7)
ggsave("cov5.png", plot = cov5, path = model_dir, width = 9, height = 7)

setwd(model_dir)
plots_cov <- list(cov2, cov3, cov4, cov5)
save(plots_cov, file = "cov_500.RData" )


se2 <- plot_se(se_mean,2, titles)
se3 <- plot_se(se_mean,3, titles)
se4 <- plot_se(se_mean,4, titles)
se5 <- plot_se(se_mean,5, titles)

setwd(here())
ggsave("se2.png", plot = se2, path = model_dir, width = 9, height = 7)
ggsave("se3.png", plot = se3, path = model_dir, width = 9, height = 7)
ggsave("se4.png", plot = se4, path = model_dir, width = 9, height = 7)
ggsave("se5.png", plot = se5, path = model_dir, width = 9, height = 7)

setwd(model_dir)
plots_se <- list(se2, se3, se4, se5)
save(plots_se, file = "se_500.RData" )
setwd(here())

distr2 <- plot_distribution(estimates,2)
distr3 <- plot_distribution(estimates,3)
distr4 <- plot_distribution(estimates,4)
distr5 <- plot_distribution(estimates,5)

setwd(model_dir)
plots_distr <- list(distr2, distr3, distr4, distr5)
save(plots_distr, file = "distr_500.RData" )
setwd(here())

# 
#plfe <- plot_lfe(0)
#plfe_dem <- plot_lfe(1) 
# 

relbias_lfe <- function(df, gt_tls) {
  df <- as.data.frame(df)
  gt_vec <- rep(gt_tls, times = nrow(df) / 2)
  dem_vec <- rep(c(0,1),times = nrow(df) / 2)
  df$dem <- dem_vec
  df$lfe <- (as.numeric(as.character(df$lfe)) - gt_vec) / gt_vec
  return(df)
}

lfe_estimates <- lapply(lfe_estimates, relbias_lfe, gt_tls = gt_tls)

process_ci <- function(temp) {
  ci_lfe <- temp %>%
    group_by(model, dem) %>%
    summarise(
      across(
        lfe,
        list(
          mean = ~ round(mean(.x, na.rm = TRUE), 5),
          lower_ci = ~ round(
            mean(.x, na.rm = TRUE) - qt(0.975, df = sum(!is.na(.x)) - 1) * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),5),
          upper_ci = ~ round(
            mean(.x, na.rm = TRUE) + qt(0.975, df = sum(!is.na(.x)) - 1) * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),5)
          )),
      .groups = 'drop'
    )
  
  desired_order <- c("flexsurv_EO", "coxph", "flexsurv", "msm", "msm_age", "nhm", "imputation")
  ci_lfe$model <- factor(ci_lfe$model, levels = desired_order, ordered = TRUE)
  ci_lfe <- ci_lfe[order(ci_lfe$model), ]
  
  return(ci_lfe)
}

ci_lfe <- lapply(lfe_estimates, process_ci)


plfe <- plot_lfe_bias(ci_lfe)

setwd(model_dir)
save(plfe, file = "plfe500.RData" )


significant_covs <- data.frame("cov1"= c(0,1,1), "cov2"= c(1,1,0), "cov3"=c(1,1,1), "transition"=c(1,2,3))

pw2 <- power_categorical(significant_covs, significancy, scheme=2)
pw3 <- power_categorical(significant_covs, significancy, scheme=3)
pw4 <- power_categorical(significant_covs, significancy, scheme=4)
pw5 <- power_categorical(significant_covs, significancy, scheme=5)

setwd(model_dir)
plots_power <- list(pw2[[1]], pw2[[2]], pw2[[3]], pw3[[1]], pw3[[2]], pw3[[3]], pw4[[1]], pw4[[2]], pw4[[3]], pw5[[1]], pw5[[2]], pw5[[3]])
save(plots_power, file = "power_500.RData" )
setwd(here())

err2 <- type_1_error(significant_covs, significancy, scheme=2)
err3 <- type_1_error(significant_covs, significancy, scheme=3)
err4 <- type_1_error(significant_covs, significancy, scheme=4)
err5 <- type_1_error(significant_covs, significancy, scheme=5)

setwd(model_dir)
plots_errorI <- list(err2,err3,err4,err5)
save(plots_errorI, file = "typeIerr_500.RData" )
setwd(here())


computational_time <- rbind(ct_all_schemes[[1]], ct_all_schemes[[2]], ct_all_schemes[[3]], ct_all_schemes[[4]])
colnames(computational_time) <- c("a","b","c","d","e","f")
computational_time <- as.data.frame(computational_time)
computational_time$scheme <- rep(2:5, each=100)
mean_ct <- computational_time %>%
  group_by(scheme)%>%
  summarise(across(everything(), mean, na.rm = TRUE))
mean_ct <- mean_ct%>%
  mutate(
    scheme = case_when(
      scheme == "2" ~ "1 year",
      scheme == "3" ~ "3 years",
      scheme == "4" ~ "3-6 years",
      scheme == "5" ~ "irregular",
      ))
setwd(model_dir)
save(mean_ct, file = "ct_2K.RData" )
setwd(here())


# titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
# plot1 <- plot_ct(2, titles, combined_cov[[1]])
# plot2 <- plot_ct(3, titles, combined_cov[[2]])
# plot3 <- plot_ct(4, titles, combined_cov[[3]])
# plot4 <- plot_ct(5, titles, combined_cov[[4]])
# 
# ct_plot <- (plot1 + theme(legend.position = "none")) +
#   (plot2 + theme(legend.position = "none")) +
#   (plot3 + theme(legend.position = "none")) +
#   (plot4 + theme(legend.position = "none")) +       
#   plot_layout(guides = "collect")  
# 
# ct_plot <- ct_plot & theme(legend.position = "left")
# 
# ggsave("ct.png", plot = ct_plot, path = model_dir, width = 9, height = 7)






