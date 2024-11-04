compute_bias_rel <- function( model_params, ground_truth_params){
  bias <- matrix(0, nrow = 3, ncol = ncol(ground_truth_params))
  colnames(bias) <- colnames(ground_truth_params)
  rownames(bias) <- rownames(ground_truth_params)
  
  common_cols <- intersect(colnames(model_params), colnames(ground_truth_params))
  common_indices_truth <- match(common_cols, colnames(ground_truth_params))
  common_indices_params <- match(common_cols, colnames(model_params))
  
  for (i in seq_along(common_indices_truth)) {
    for (j in 1:3){
      if (ground_truth_params[j, common_indices_truth[i]] != 0) {
        bias[j, common_indices_truth[i]] <- (model_params[j, common_indices_params[i]] - ground_truth_params[j, common_indices_truth[i]]) / ground_truth_params[j, common_indices_truth[i]]
      } else {
        bias[j, common_indices_truth[i]] <- NA
      }
    }
  }
  
  return(bias)
}