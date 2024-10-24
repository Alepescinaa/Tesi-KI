get_params_nhm <- function (x, ci = TRUE, ...) 
{
  model <- x
  par <- model$par
  fisher <- model$hess
  par_names <- model$parnames
  if (is.null(fisher)) 
    fisher <- model$fisher
  if (is.null(fisher) & ci) 
    stop("Fisher information required for confidence intervals")
  if (is.null(par_names)) 
    par_names <- sapply(1:length(par), function(x) paste("Parameter", 
                                                         x))
  if (model$singular & ci) {
    warning("Cannot calculate confidence intervals due to a singular Hessian")
    ci <- FALSE
  }
  if (ci) {
    std.err <- diag(solve(fisher))^0.5
    mat <- round(cbind(par, par - qnorm(0.975) * std.err, 
                       par + qnorm(0.975) * std.err), 4)
    dimnames(mat)[[2]] <- c("Est", "Low 95%", "Up 95%")
  }
  else {
    mat <- round(cbind(par), 4)
    dimnames(mat)[[2]] <- c("Est")
  }
  if (!is.null(model$fixedpar)) {
    mat2 <- array(NA, c(model$npar, dim(mat)[2]))
    mat2[-model$fixedpar, ] <- mat
    mat2[model$fixedpar, 1] <- model$fixval
    dimnames(mat2)[[2]] <- dimnames(mat)[[2]]
    mat <- mat2
  }
  dimnames(mat)[[1]] <- par_names
  return(mat)
  
}
