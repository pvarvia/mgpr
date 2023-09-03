# A function for kfold cross validation
kfold_mgpr <- function(mgpr, newdatax = NULL, credinter = NULL,
                       covout = FALSE, kfold = 0, fixneg = FALSE,
                       verbose = FALSE, seed = NULL) {
  # Check arguments
  if (missing(mgpr) | !is.list(mgpr)) {
    stop("Error: mgpr is not a list or missing.")
  }
  
  if (!is.logical(verbose)) {
    stop("Invalid verbose argument.", fill = TRUE)
  }
  
  if (!is.logical(covout)) {
    stop("Error: covout must be boolean.")
  }
  
  if (!is.logical(fixneg)) {
    stop("Error: Invalid fixneg argument. Must be boolean.")
  }
  
  if (missing(kfold) | kfold == 0 || !is.numeric(kfold)) {
    stop(paste0(
      "Please set group size (k) parameter! When k = nrow(data) ",
      "= Leave-one-out cross validation"
    ))
  } else {
    k <- kfold
    kfold <- NULL
    if (k == 1) {
      stop("Error: k = 1 is not allowed. Use kfold = 0 or kfold > 1.")
    }
    if (k == dim(mgpr$trainy)[1]) {
      if (verbose) {
        cat("leave-one-out cross validation", fill = TRUE)
      }
      kfold <- FALSE
    } else if (dim(mgpr$trainy)[1] %% k == 0) {
      if (verbose) {
        cat("k-fold cross validation with k = ", k, fill = TRUE)
      }
      kfold <- TRUE
    } else {
      if (verbose) {
        cat(paste0(
          "WARNING: nrow(data) %% k != 0. The folds are not ",
          "equal-sized."
        ), fill = TRUE)
      }
      kfold <- TRUE
    }
  }
  
  if (!is.null(credinter) && !is.numeric(credinter)) {
    stop("Error: Invalid credinter argument.")
  } else {
    if (is.numeric(credinter)) {
      if (credinter <= 0 | credinter >= 1) {
        stop(paste0(
          "Error: Define credible intervals! ",
          "An example: set 0.95 to get 95% credible intervals"
        ))
      }
    }
  }
  
  if (!is.null(newdatax)) {
    stop("Error: newdata is not allowed because of kfold.")
  }
  #####################################
  error_variance <- mgpr$errorvar
  datax <- t(mgpr$trainx)
  datay <- t(mgpr$trainy)
  n_n <- dim(datay)[2] # predictor n
  ny <- dim(mgpr$trainy)[2] # response n
  # predicted attributes go here
  xgp <- array(rep(0, ny * n_n), dim = c(ny, n_n))
  
  # Initialize the following matrices if CIs are computed
  if (!is.null(credinter)) {
    # Credible intervals
    lowerlimgp <- array(rep(0, ny * n_n), dim = c(ny, n_n))
    upperlimgp <- array(rep(0, ny * n_n), dim = c(ny, n_n))
  }
  # If covariance output requested
  if (covout) {
    # Covariances
    gpcovar <- array(rep(0, ny * ny * n_n), dim = c(ny, ny, n_n))
  }
  
  # set seed if given
  if (!is.null(seed)) {
    oldstate <- .Random.seed
    set.seed(seed)
  }
  
  #####################################
  foldcount_user <- k
  
  if (kfold) {
    kfold_sel <- seq(1, n_n, 1)
    # Remainder will be added
    remainder <- n_n %% k
    foldsize_user <- (n_n - remainder) / k
    remainder_per_iter <- (remainder -
                             (remainder %% foldcount_user)) / foldcount_user
    # distribute the rest randomly (remainder)
    excess_count <- remainder - (remainder_per_iter * foldcount_user)
    # decide folds
    add_excess <- sample(seq(1, foldcount_user, 1),
                         size = excess_count, replace = FALSE
    )
  }
  
  for (chosengroup in 1:foldcount_user) {
    start_time <- Sys.time()
    
    if (kfold) {
      # Fix the fold size (equally distributed)
      k <- foldsize_user + remainder_per_iter
      # Logical check; will be checked later on.
      remainder <- remainder - remainder_per_iter
      if (chosengroup %in% add_excess) {
        k <- k + 1
        remainder <- remainder - 1
      }
      if (verbose) {
        cat("k-fold cross validation: ", chosengroup, "/", foldcount_user,
            fill = TRUE
        )
        cat("Fold number: ", chosengroup, "\nFold size: ", k, fill = TRUE)
      }
      chosenplotrows <- sample(kfold_sel, k, replace = FALSE)
      # update kfold selection, i.e. remove selected rows.
      kfold_sel <- kfold_sel[-which(kfold_sel %in% chosenplotrows)]
    } else {
      # LOOCV
      chosenplotrows <- chosengroup
      if (verbose) {
        if (chosengroup %% 50 == 0) {
          cat("leave-one-out cross validation: ", chosenplotrows, "/", n_n,
              fill = TRUE
          )
        }
      }
    }
    
    # select target data
    z <- datax[, chosenplotrows]
    # remove row
    datax_lo <- datax[, -chosenplotrows]
    datay_lo <- datay[, -chosenplotrows]
    # check if only one predictor
    if (dim(mgpr$trainx)[2] == 1) {
      datax_lo <- t(array(datax_lo,
                          dim = c(dim(datax)[2] - length(chosenplotrows), 1)
      ))
      # then fix also z
      z <- (array(z, dim = c(1, length(z))))
    }
    # check if only one response
    if (dim(mgpr$Cy)[1] == 1 & dim(mgpr$Cy)[2] == 1) {
      datay_lo <- t(array(datay_lo, dim = c(length(datay_lo), 1)))
    }
    
    # remove target data rows and columns from kernel matrices
    K_lo <- mgpr$K[-chosenplotrows, -chosenplotrows]
    Ie_lo <- mgpr$Ie[-chosenplotrows, -chosenplotrows]
    
    # Update mgpr object
    mod1 <- mgpr
    mod1$K <- K_lo
    mod1$Ie <- Ie_lo
    mod1$trainy <- t(datay_lo)
    mod1$trainx <- t(datax_lo)
    mod1$trainMeanSd <- NULL
    nx1 <-  dim(t(mod1$trainy))[1]
    # Compute prior expectancy for training data
    mu_Y <- apply(t(mod1$trainx), 2, FUN = mod1$meanf)
    # precomputed parts of prediction equations
    if (nx1 == 1) {
      mod1$predM1 <- chol2inv(chol(as.vector(mod1$Cy) * mod1$K +
                                     mod1$E * mod1$Ie))
      mod1$predM2 <- mod1$predM1 %*% as.vector(mod1$trainy - mu_Y)
    } else {
      mod1$predM1 <- chol2inv(chol(kronecker(mod1$Cy, mod1$K) +
                                     kronecker(mod1$E, mod1$Ie)))
      mod1$predM2 <- mod1$predM1 %*% as.vector(mod1$trainy - t(mu_Y))
    }
    
    
    # Prediction
    if (is.null(credinter)) {
      pred <- predict(
        mgpr = mod1,
        newdatax = t(z),
        verbose = verbose,
        fixneg = fixneg,
        covout = covout
      )
      # Append predictions to correct places
      if (!is.data.frame(pred)) { # Ret a list if credinter/covout requested
        for (kfoldrow in 1:length(chosenplotrows)) {
          xgp[, chosenplotrows[kfoldrow]] <- as.numeric(pred$pred[kfoldrow, ])
          if (covout) {
            gpcovar[, , chosenplotrows[kfoldrow]] <- pred$covpred[, , kfoldrow]
          }
        }
      } else { # This case with covout = FALSE and credinter = NULL
        for (kfoldrow in 1:length(chosenplotrows)) {
          xgp[, chosenplotrows[kfoldrow]] <- as.numeric(pred[kfoldrow, ])
        }
      }
    } else {
      pred <- predict(
        mgpr = mod1,
        newdatax = t(z),
        credinter = credinter,
        verbose = verbose,
        fixneg = fixneg,
        covout = covout
      )
      
      # append predictions and CIs to correct places..
      # predict returns a list
      if (kfold) {
        for (kfoldrow in 1:length(chosenplotrows)) {
          xgp[, chosenplotrows[kfoldrow]] <- as.numeric(pred$pred[
            kfoldrow,
            1:dim(xgp)[1]
          ])
        }
        for (kfoldrow in 1:length(chosenplotrows)) {
          if (covout) {
            gpcovar[, , chosenplotrows[kfoldrow]] <- pred$covpred[, , kfoldrow]
          }
          lowerlimgp[, chosenplotrows[kfoldrow]] <-
            as.numeric(pred$credinter[kfoldrow, ])[1:dim(mod1$trainy)[2]]
          upperlimgp[, chosenplotrows[kfoldrow]] <-
            as.numeric(pred$credinter[kfoldrow, ])[(dim(mod1$trainy)[2] +
                                                      1):(dim(mod1$trainy)[2] * 2)]
        }
      } else {
        if (covout) {
          gpcovar[, , chosenplotrows] <- pred$covpred[, , 1]
        }
        lowerlimgp[, chosenplotrows] <- as.numeric(pred$credinter)[
          1:dim(mod1$trainy)[2]
        ]
        upperlimgp[, chosenplotrows] <- as.numeric(pred$credinter)[
          (dim(mod1$trainy)[2] + 1):(dim(mod1$trainy)[2] * 2)
        ]
        xgp[, chosenplotrows] <- as.numeric(pred$pred)
      }
    }
    end_time <- Sys.time()
    # if (verbose) {
    #   cat(paste0("Run time for kfold iteration: ",
    #    round(difftime(end_time,as.POSIXct(start_time),
    #    units="secs"),2),"secs \n"))
    # }
  }
  
  if (kfold) {
    if (remainder != 0) {
      stop("Error: k-fold cv.", remainder)
    }
  }
  
  # restore RNG state if seed was set
  if (!is.null(seed)) {
    .Random.seed <- oldstate
  }
  
  # cast to data frame and transpose to output (n x n_Y)
  xgp <- as.data.frame(t(xgp))
  if (mgpr$YColName[1] != -9999) {
    colnames(xgp) <- mgpr$YColName
  }
  if (mgpr$XRowName[1] != -9999) {
    rownames(xgp) <- mgpr$XRowName
  }
  
  # Pred covpreds for output
  if (covout) {
    # set dim names
    gpcovar <- provideDimnames(gpcovar,
                               base = list(
                                 colnames(xgp),
                                 colnames(xgp), paste(rep(
                                   "plot_",
                                   dim(gpcovar)[3]
                                 ), seq(1, dim(gpcovar)[3], 1),
                                 sep = ""
                                 )
                               )
    )
  }
  # Returns
  if (is.null(credinter)) {
    if (covout) {
      out <- list(pred = xgp, covpred = gpcovar)
    } else {
      out <- as.data.frame(xgp)
    }
    return(out)
  } else {
    CredInter <- data.frame(t(lowerlimgp), t(upperlimgp))
    colnames(CredInter) <- c(
      paste(colnames(xgp), "_lowerlimCI", sep = ""),
      paste(colnames(xgp), "_upperlimCI", sep = "")
    )
    out <- list(pred = xgp, credinter = CredInter)
    if (covout) {
      out[["covpred"]] <- gpcovar
    }
    return(out)
  }
}
