####################################
# Upload library and data
####################################

library(fs)

main_dir <- "/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/wrapper_MM/saved_models_scheme"
scheme <-  2
seed <- 1

if (scheme == 2){
  scheme_dir <- "/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/wrapper_MM/saved_models_scheme2"
} else if (scheme == 3){
  scheme_dir <- "/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/wrapper_MM/saved_models_scheme3"
} else if (scheme == 4){
  scheme_dir <- "/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/wrapper_MM/saved_models_scheme4"
} else if (scheme == 5){
  scheme_dir <- "/Users/AlessandraPescina/OneDrive - Politecnico di Milano/ANNO 5/secondo semestre/TESI/Tesi/Tesi_code/wrapper_MM/saved_models_scheme5"
}
  

seed_dir <- file.path(scheme_dir, paste0("seed_", seed))

if (dir.exists(seed_dir)) {
  setwd(seed_dir)
  
  files_to_load <- c("cox_model.RData", 
                     "flexsurv_model.RData", 
                     "msm_model.RData", 
                     "model_msm_age.RData", 
                     "model_nhm.RData", 
                     "results_imp.RData", 
                     "computational_time.RData")
  for (file in files_to_load) {
    if (file.exists(file)) {
      load(file)
    } else {
      warning(paste("File does not exist:", file))
    }
  }
} else {
  warning(paste("Seed directory does not exist:", seed_dir))
}

