compute_bias <- function( model_params, ground_truth_params){
  bias <- matrix(0, nrow = 3, ncol = ncol(ground_truth_params))
  colnames(bias) <- colnames(ground_truth_params)
  rownames(bias) <- rownames(ground_truth_params)
  
  common_cols <- intersect(colnames(model_params), colnames(ground_truth_params))
  common_indices_truth <- match(common_cols, colnames(ground_truth_params))
  common_indices_params <- match(common_cols, colnames(model_params))
  
  
  bias[, common_indices_truth] <- abs(model_params[, common_indices_params] - ground_truth_params[, common_indices_truth])
  
  return(bias)
}