totlos.simfs.mine <- function (x, trans, t = 1, tstart = 0, start = 1, newdata = NULL, ci = FALSE, 
                               tvar = "trans", tcovs = NULL, group = NULL, M = 1e+05, B = 1000, 
                               cl = 0.95, cores = NULL) 
{
  if (length(t) > 1) 
    stop("\"t\" must be a single number")
  
  # Simulate data
  sim <- sim.fmsm.mine(x = x, trans = trans, t = t, newdata = newdata, 
                  start = start, tstart= t_start, M = M, tvar = tvar, tcovs = tcovs, debug = FALSE)
  

  # Calculate the time differences and total length of stay
  dt <- diff(t(cbind(sim$t, t)))
  st <- factor(t(sim$st), levels = 1:nrow(trans))
  res <- tapply(dt, st, sum) / M
  res[is.na(res)] <- 0
  
  # Grouping states if 'group' is specified
  if (!is.null(group)) {
    if (length(group) != nrow(trans)) 
      stop("\"group\" must be a vector of length ", nrow(trans), " = number of states")
    res <- tapply(res, group, sum)
  }
  
  # Confidence intervals if 'ci' is TRUE
  if (ci) {
    resci <- bootci.fmsm(x, B, fn = totlos.simfs, ci = FALSE, 
                         cl = cl, cores = cores, trans = trans, t = t, 
                         start = start, newdata = newdata, tvar = tvar, 
                         tcovs = tcovs, group = group, M = M, tstart = tstart)
    resl <- resci[1, ]
    resu <- resci[2, ]
    names(resl) <- names(resu) <- t
    attr(res, "lower") <- resl
    attr(res, "upper") <- resu
    class(res) <- "fs.msm.est"
    res <- cbind(est = res, L = resl, U = resu)
  }
  
  res
}
