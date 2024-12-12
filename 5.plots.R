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
  "./wrapper_MM/functions_performance/plot_lfe_bias.R"
  
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

ggsave("convergence.png", plot = combined_plot, path = model_dir, width = 10, height = 7) 


titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
pb2 <- plot_bias(2, titles)
pb3 <- plot_bias(3, titles)
pb4 <- plot_bias(4, titles)
pb5 <- plot_bias(5, titles)

ggsave("bias2.png", plot = pb2, path = model_dir, width = 9, height = 7)
ggsave("bias3.png", plot = pb3, path = model_dir, width = 9, height = 7)
ggsave("bias4.png", plot = pb4, path = model_dir, width = 9, height = 7)
ggsave("bias5.png", plot = pb5, path = model_dir, width = 9, height = 7)


titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
cov2 <- plot_coverage(2, titles)
cov3 <- plot_coverage(3, titles)
cov4 <- plot_coverage(4, titles)
cov5 <- plot_coverage(5, titles)

ggsave("cov2.png", plot = cov2, path = model_dir, width = 9, height = 7)
ggsave("cov3.png", plot = cov3, path = model_dir, width = 9, height = 7)
ggsave("cov4.png", plot = cov4, path = model_dir, width = 9, height = 7)
ggsave("cov5.png", plot = cov5, path = model_dir, width = 9, height = 7)


#da sistemare per rappresentare dev standard
# titles <- c("Population Based Study (1 year)", "Population Based Study (3 years)", "Population Based Study (3-6 years)", "Electronic Health Record")
# sd2 <- plot_width(mean_width,2, titles)
# sd3 <- plot_width(mean_width,3, titles)
# sd4 <- plot_width(mean_width,4, titles)
# sd5 <- plot_width(mean_width,5, titles)
# 
# ggsave("sd2.png", plot = sd2, path = model_dir, width = 9, height = 7)
# ggsave("sd3.png", plot = sd3, path = model_dir, width = 9, height = 7)
# ggsave("sd4.png", plot = sd4, path = model_dir, width = 9, height = 7)
# ggsave("sd5.png", plot = sd5, path = model_dir, width = 9, height = 7)


# add new version of error type one and power

plfe <- plot_lfe(0)
plfe_dem <- plot_lfe(1) 

plfe_b <- plot_lfe_bias(0)
plfe_dem_b <- plot_lfe_bias(1)

ggsave("lfe.png", plot = plfe, path = model_dir, width = 10, height = 7) 
ggsave("years_dem.png", plot = plfe_dem, path = model_dir, width = 10, height = 7) 
ggsave("lfe_b.png", plot = plfe_b, path = model_dir, width = 10, height = 7) 
ggsave("years_dem_b.png", plot = plfe_dem_b, path = model_dir, width = 10, height = 7) 


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

ggsave("ct.png", plot = ct_plot, path = model_dir, width = 9, height = 7)





