is.flexsurvlist <- function (x) {
  is.list(x) && (length(x) > 0) && inherits(x[[1]], "flexsurvreg") &&
    all(sapply(x, inherits, "flexsurvreg"))
}
