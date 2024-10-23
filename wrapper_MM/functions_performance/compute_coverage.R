compute_coverage <- function(lower_bounds, upper_bounds, ground_truth_params) {
  coverage <- matrix(FALSE, nrow = nrow(ground_truth_params), ncol = ncol(ground_truth_params))
  colnames(coverage) <- colnames(ground_truth_params)
  rownames(coverage) <- rownames(ground_truth_params)
  
  common_cols <- intersect(rownames(lower_bounds), colnames(ground_truth_params))
  
  for (col in common_cols) {
    truth_index <- which(colnames(ground_truth_params) == col)
    lower_index <- which(rownames(lower_bounds) == col)
    upper_index <- which(rownames(upper_bounds) == col)
    
    coverage[, truth_index] <- as.numeric(
      ground_truth_params[, truth_index] >= lower_bounds[lower_index, ] &
        ground_truth_params[, truth_index] <= upper_bounds[upper_index, ]
    )
  }
  
  return(coverage)
}
