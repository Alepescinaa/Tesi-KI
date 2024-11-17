smms_mine <- function (startval, data, graph, X = NULL, abs_exact = TRUE, 
          mc_cores = 1, hessian_matrix = FALSE, cmethod = "hcubature") 
{
  formula_obs_types = all_types(graph)
  edge_mats = edge_matrices(graph)
  state_ord = state_ordering(graph)
  names_surv_dens = names_of_survival_density(graph)
  absorbing_states <- sort(state_ord$order[which(state_ord$type == 
                                                   "abs")])
  all_data_set = arrange_data(data, graph)
  timepointMat <- all_data_set[, 1:(dim(all_data_set)[2] - 
                                      1)]
  observation_type = rep(NA, nrow(all_data_set))
  all_integral_limits = list()
  integrand = list()
  for (i in 1:nrow(all_data_set)) {
    observation_type[i] = all_data_set[i, "obs_type"]
    f_types = names(which(formula_obs_types[, observation_type[i]] == 
                            1))
    integrand_mellomregn = list()
    integral_mellomregn = list()
    for (j in 1:length(f_types)) {
      integrand_mellomregn[[j]] = eval(parse(text = type_to_integrand(f_types[j], 
                                                                      edge_mats, names_surv_dens, abs_exact = abs_exact)))
      integral_mellomregn[[j]] = finding_limits(timepointMat[i, 
      ], f_types[j], edge_mats, absorbing_states, 
      abs_exact = abs_exact)
    }
    all_integral_limits[[i]] = integral_mellomregn
    integrand[[i]] = integrand_mellomregn
  }
  optimizer <- stats::optim(startval,mloglikelihood,integrand = integrand,limits = all_integral_limits,X=X, method = "L-BFGS",
                            lower = rep(-50,length(startval)), upper = rep(50, length(startval)), mc_cores=mc_cores,hessian = FALSE)
  
  if (hessian_matrix == TRUE) {
    hessian_optimizer = numDeriv::hessian(mloglikelihood, 
                                          optimizer$par, integrands = integrand, limits = all_integral_limits, 
                                          cmethod = cmethod, mc_cores = mc_cores, X = X)
    return(list(opt = optimizer, hess = hessian_optimizer))
  }
  else {
    return(list(opt = optimizer))
  }
}