init_covf <- function(covftype) {
  if (identical(covftype, "matern12")) {
    # Matern 1/2 a.k.a Ornstein-Uhlenbeck
    # d = distance
    covfunc <- function(d, ksigma, corlen) {
      ksigma^2 * exp(-d / corlen)
    }
  } else if (identical(covftype, "matern32")) {
    # Matern 3/2
    covfunc <- function(d, ksigma, corlen) {
      ksigma^2 * (1 + sqrt(3) *
                    d / corlen) * exp(-sqrt(3) * d / corlen)
    }
  } else if (identical(covftype, "matern52")) {
    # Matern 5/2
    covfunc <- function(d, ksigma, corlen) {
      ksigma^2 * (1 + sqrt(5) *
                    d / corlen + 5 * (d^2) / (3 * corlen^2)) *
        exp(-sqrt(5) * d / corlen)
    }
  } else if (identical(covftype, "rbf")) {
    # RBF, Gaussian kernel, Squared exponential
    covfunc <- function(d, ksigma, corlen) {
      ksigma^2 * exp(-(d^2) / (2 * corlen^2))
    }
  } else {
    stop("User-defined kernel is not supported in this version.")
  }
  
  
}