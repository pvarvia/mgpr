init_meanf <- function(meanftype,mgpr,verbose) {
  
  if (!(identical(meanftype, "zero") ||
        identical(meanftype, "avg") ||
        identical(meanftype, "linear"))) {
    stop(paste0(
      "Error: Invalid meanf argument. ",
      "Type 'zero', 'avg', or 'linear'"
    ))
  }
  
  if (!is.function(meanftype)) {
    nx <- dim(mgpr$trainy)[2]
    
    # Zero mean
    if (identical(meanftype, "zero")) {
      #remove unneeded mgpr-object to avoid duplicate
      rm(mgpr)
      meanf <- function(x) {
        # Return
        return(rep(0, nx))
      }
    } else if (identical(meanftype, "avg")) { # sample average of training data
      mux <- t(apply(mgpr$trainy, MARGIN = 2, FUN = mean))
      rm(mgpr)
      meanf <- function(x) {
        # Return
        return(mux)
      }
    } else if (identical(meanftype, "linear")) { # linear estimates
      # No need for special cases for single response case since
      # the format is always ok in mgpr object
      mux <- t(apply(mgpr$trainy, MARGIN = 2, FUN = mean))
      Ccmf <- cov(cbind(mgpr$trainx, mgpr$trainy))
      muz <- t(apply(mgpr$trainx, MARGIN = 2, FUN = mean))
      Cxz <- t(Ccmf[
        1:(nrow(Ccmf) - (nx)),
        (ncol(Ccmf) - (nx - 1)):ncol(Ccmf)
        ])
      Cz <- Ccmf[1:(nrow(Ccmf) - nx), 1:(ncol(Ccmf) - nx)]
      # Cz should be at least psd, solve using Cholesky
      Acm <- try(Cxz %*% chol2inv(chol(Cz)), silent = TRUE)
      if (identical(class(Acm), "try-error")) {
        if (verbose) {
          cat(paste0(
            "Warning: Do you have highly correlated features? ",
            "System is computationally singular. Using",
            " regularized inverse of Cz matrix.",
            "Using the linear mean function."
          ),
          fill = TRUE
          )
        }
        Acm <- try(Cxz %*% chol2inv(chol(Cz + 1e-8 * diag(dim(mgpr$trainx)[2]))))
      }
      rm(mgpr)
      meanf <- function(x) {
        # Return
        return(mux + t(Acm %*% array(x - muz, dim = c(length(x), 1))))
      }
    }
  } else {
    stop("Error: User-defined mean function is not supported in this version.")
  }
  
  
}