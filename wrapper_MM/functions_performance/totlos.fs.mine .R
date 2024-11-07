totlos.fs.mine <- function (x, trans = NULL, t_start= 0, t = 1, newdata = NULL, ci = FALSE, 
                            tvar = "trans", sing.inf = 1e+10, B = 1000, cl = 0.95, ...) 
{
  if (is.null(trans)) {
    if (!is.null(attr(x, "trans"))) 
      trans <- attr(x, "trans")
    else stop("`trans` not supplied and not found in `x`")
  }
  ntr <- sum(!is.na(trans))
  n <- nrow(trans)
  nsq <- n * n
  dp <- function(t, y, parms, ...) {
    P <- matrix(y[nsq + 1:nsq], nrow = n, ncol = n)
    haz <- numeric(n)
    for (i in 1:ntr) {
      xi <- if (is.flexsurvlist(x)) 
        x[[i]]
      else x
      hcall <- list(x = t)
      for (j in seq_along(xi$dlist$pars)) hcall[[xi$dlist$pars[j]]] <- parms$par[[i]][j]
      hcall <- c(hcall, parms$aux[[i]])
      haz[i] <- do.call(xi$dfns$h, hcall)
    }
    Q <- haz[trans]
    Q[is.na(Q)] <- 0
    Q[is.infinite(Q) & Q > 0] <- sing.inf
    Q <- matrix(Q, nrow = n, ncol = n)
    diag(Q) <- -rowSums(Q)
    list(cbind(P, P %*% Q))
  }
  nt <- length(t)
  if (nt < 1) 
    stop("number of times should be at least one")
  basepar <- flexsurv::pars.fmsm(x = x, trans = trans, newdata = newdata, 
                       tvar = tvar)
  auxpar <- if (!is.flexsurvlist(x)) 
    x$aux
  else lapply(x, function(x) x$aux)
  init <- cbind(matrix(0, nrow = n, ncol = n), diag(n))
  res <- ode(y = init, times = c(t_start, t), func = dp, parms = list(par = basepar, 
                                                                aux = auxpar), ...)[-1, -1]
  res.t <- lapply(split(res, 1:nt), function(x) matrix(x[1:nsq], 
                                                       nrow = n))
  res.p <- lapply(split(res, 1:nt), function(x) matrix(x[nsq + 
                                                           1:nsq], nrow = n))
  names(res.t) <- names(res.p) <- t
  if (ci) {
    resci <- flexsurv::bootci.fmsm(x, B, fn = totlos.fs, attrs = "P", 
                         cl = cl, ci = FALSE, trans = trans, t = t, newdata = newdata, 
                         tvar = tvar, sing.inf = sing.inf, ...)
    tind <- rep(rep(1:nt, each = n * n), 2)
    res.tl <- lapply(split(resci[1, ], tind), function(x) matrix(x[1:nsq], 
                                                                 nrow = n))
    res.tu <- lapply(split(resci[2, ], tind), function(x) matrix(x[1:nsq], 
                                                                 nrow = n))
    res.pl <- lapply(split(resci[1, ], tind), function(x) matrix(x[nsq + 
                                                                     1:nsq], nrow = n))
    res.pu <- lapply(split(resci[2, ], tind), function(x) matrix(x[nsq + 
                                                                     1:nsq], nrow = n))
    names(res.tl) <- names(res.tu) <- names(res.pl) <- names(res.pu) <- t
    for (i in 1:nt) {
      attr(res.t[[i]], "lower") <- res.tl[[i]]
      attr(res.t[[i]], "upper") <- res.tu[[i]]
      class(res.t[[i]]) <- "fs.msm.est"
      attr(res.p[[i]], "lower") <- res.pl[[i]]
      attr(res.p[[i]], "upper") <- res.pu[[i]]
      class(res.p[[i]]) <- "fs.msm.est"
    }
  }
  if (nt == 1) {
    res.t <- res.t[[1]]
    res.p <- res.p[[1]]
  }
  attr(res.t, "P") <- res.p
  class(res.t) <- "totlos.fs"
  res.t
}