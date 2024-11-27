transient <- function(trans){
  which(apply(trans, 1, function(x)any(!is.na(x))))
}

is.flexsurvlist <- function (x) {
  is.list(x) && (length(x) > 0) && inherits(x[[1]], "flexsurvreg") &&
    all(sapply(x, inherits, "flexsurvreg"))
}

sim.fmsm.mine <- function (x, trans = NULL, t, newdata = NULL, start = 1, tstart=t_start, M = 10, 
  tvar = "trans", tcovs = NULL, tidy = FALSE, debug = FALSE) 
{
  if (is.null(trans)) {
    if (!is.null(attr(x, "trans"))) 
      trans <- attr(x, "trans")
    else stop("`trans` not supplied and not found in `x`")
  }
  if (length(t) == 1) 
    t <- rep(t, M)
  else if (length(t) != M) 
    stop("length of t should be 1 or M=", M)
  if (length(start) == 1) 
    start <- rep(start, M)
  else if (length(start) != M) 
    stop("length of start should be 1 or M=", M)
  basepars.all <- pars.fmsm(x = x, trans = trans, newdata = newdata, 
    tvar = tvar)
  nst <- nrow(trans)
  res.st <- cur.st <- start
  res.t <- cur.t <- rep(0, M)
  todo <- seq_len(M)
  while (any(todo)) {
    if (debug) {
      cat("cur.t\n")
      cat(cur.t)
      cat("\n")
    }
    if (debug) {
      cat("cur.st\n")
      cat(cur.st)
      cat("\n")
    }
    if (debug) {
      cat("TODO\n")
      cat(todo)
      cat("\n")
    }
    cur.st.out <- cur.st[todo]
    cur.t.out <- cur.t[todo]
    done <- numeric()
    for (i in unique(cur.st[todo])) {
      if (i %in% transient(trans)) {
        transi <- na.omit(trans[i, ])
        ni <- sum(cur.st[todo] == i)
        t.trans1 <- matrix(0, nrow = ni, ncol = length(transi))
        for (j in seq_along(transi)) {
          xbase <- if (is.flexsurvlist(x)) 
            x[[transi[j]]]
          else x
          if (length(tcovs) > 0) {
            basepars <- form.basepars.tcovs(x, transi[j], 
              newdata, tcovs, cur.t.out)
          }
          else basepars <- as.list(as.data.frame(basepars.all[[transi[j]]]))
          fncall <- c(list(n = ni), basepars, xbase$aux)
          if (is.null(xbase$dfns$r)) 
            stop("No random sampling function found for this model")
          t.trans1[, j] <- do.call(xbase$dfns$r, fncall)
        }
        if (debug) {
          print(t(t.trans1))
        }
        mc <- max.col(-t.trans1)
        if (debug) {
          cat("mc\n")
          cat(mc)
          cat("\n")
        }
        if (debug) {
          cat("transi\n")
          cat(transi)
          cat("\n")
        }
        if (debug) {
          cat("trans[i,]\n")
          cat(trans[i, ])
          cat("\n")
        }
        next.state <- match(transi[mc], trans[i, ])
        if (debug) {
          cat("next.state\n")
          cat(next.state)
          cat("\n")
        }
        next.time <- t.trans1[cbind(seq_along(next.state), 
          mc)]
        if (debug) {
          cat("next.time\n")
          cat(next.time)
          cat("\n")
        }
        inds <- which(cur.st[todo] == i)
        cur.t.out[inds] <- cur.t.out[inds] + next.time
        cens <- cur.t.out[inds] > t[inds]
        cur.t.out[inds][cens] <- t[inds][cens]
        cur.st.out[inds][!cens] <- next.state[!cens]
        done <- todo[inds][cens]
      }
    }
    cur.st[todo] <- cur.st.out
    cur.t[todo] <- cur.t.out
    res.st <- cbind(res.st, cur.st)
    res.t <- cbind(res.t, cur.t)
    done <- union(done, which(cur.st %in% absorbing(trans)))
    todo <- setdiff(todo, done)
    if (debug) {
      cat("\n")
    }
  }
  res <- list(st = unname(res.st), t = unname(res.t))
  attr(res, "trans") <- trans
  attr(res, "statenames") <- attr(x, "statenames")
  if (tidy) 
    res <- simfs_bytrans(res)
  res
}