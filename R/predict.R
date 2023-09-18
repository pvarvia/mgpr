#' Predict function for an mgpr model
#'
#' Predict using a Multivariate Gaussian process regression (mgpr) model. The
#' function supports k-fold cross validation, new predictor variables, limiting
#' predictions to positive scale, and the construction of credible intervals.
#' @param mgpr a \code{mgpr} model object, usually, a result of a call to
#' \code{\link{mgpr}}.
#' @param newdatax a data frame or an object coercible by as.data.frame to a
#' data frame containing new predictor variables. Column names must match with
#' the ones in the \code{mgpr} object. If NULL \code{datax} used in model
#' training is used.
#' @param credinter threshold of the credible intervals, 
#' e.g. \code{credinter=0.95}.
#' If NULL credible intervals are not returned.
#' @param covout logical. Determines whether the covariance matrix is returned
#' or not.
#' @param kfold value >=2 enables the k-fold cross validation mode. The
#' \code{kfold} argument defines the k value. The value must be in the range
#' 2 and nrow(data).
#' @param fixneg logical. Determines whether negative predictions (and credible
#' intervals) are corrected to be positive or zero, see details in Varvia et al.
#' (2019).
#' @param verbose logical. Report extra information on progress.
#' @return Returns predicted values. If credible intervals/covariance matrix
#' is/are requested, the output will be a list.
#' @references Varvia, P., Lahivaara, T., Maltamo, M., Packalen, P. and Seppanen,
#' A. 2019. Gaussian process regression for forest attribute estimation from
#' airborne laser scanning data. 
#' IEEE Transactions on Geoscience and Remote Sensing 57(6): 3361-3369.
#'
#' R package article to appear.
#'
#' @examples
#' data(mgprdata)
#' m <- mgpr(
#'   datay = mgprdata[, 1:3], datax = mgprdata[, 4:39],
#'   kernel = "matern32", kernpar = list(sigma = 1, corlen = 5, errorvar = 0.1)
#' )
#' p <- predict(m, credinter = 0.95)
#' str(p)
#' @export
predict <- function(mgpr,
                    newdatax = NULL,
                    credinter = NULL,
                    covout = FALSE,
                    kfold = 0,
                    fixneg = FALSE,
                    verbose = FALSE) {
  if (missing(mgpr) || !is.list(mgpr)) {
    stop("Error: mgpr is not a list or missing.")
  }
  
  if (kfold == 0) {
    UseMethod("predict")
  } else {
    if (!is.null(newdatax)) {
      stop("Error: set newdatax = NULL or kfold = 0.")
    }
    kfold_mgpr(
      mgpr = mgpr,
      newdatax = newdatax,
      credinter = credinter,
      covout = covout,
      kfold = kfold,
      fixneg = fixneg,
      verbose = verbose
    )
  }
}

#' @method predict mgpr
#' @export
predict.mgpr <- function(mgpr,
                         newdatax = NULL,
                         credinter = NULL,
                         covout = FALSE,
                         kfold = 0,
                         fixneg = FALSE,
                         verbose = FALSE) {
  Traindatax <- NULL
  Targetdatax <- NULL
  ##############################################
  # Define and check parameters
  if (missing(mgpr) || !is.list(mgpr)) {
    stop("Error: mgpr is not a list or missing.")
  }
  
  if (!is.logical(verbose)) {
    stop("Invalid verbose argument.", fill = TRUE)
  }
  
  if (!is.logical(covout)) {
    stop("Error: covout must be boolean.")
  }
  
  if (!is.numeric(kfold) || (kfold != 0)) {
    stop("Invalid kfold argument.", fill = TRUE)
  }
  
  if (!is.logical(fixneg)) {
    stop("Error: Invalid fixneg argument. Must be boolean.")
  }
  
  if (is.null(newdatax)) {
    Traindatax <- t(mgpr$trainx)
    Targetdatax <- t(mgpr$trainx)
    newdatax.rownames <- -9999
  } else {
    # Get trainx from model object
    Traindatax <- t(mgpr$trainx)
    # Collect newdata rownames if data.frame
    if (is.data.frame(newdatax)) {
      newdatax.rownames <- row.names(newdatax)
    } else {
      newdatax.rownames <- -9999
    }
    # Select predictors by column names (if data frame is used,
    # the target predictors should match with training predictors)
    # If user-input is matrix; check cannot be carried out.
    if (is.data.frame(newdatax) && mgpr$XColName[1] != -9999) {
      if (any(!(colnames(newdatax) %in% mgpr$XColName))) {
        stop(paste0(
          "Error. Predictor names in new data differ from those ",
          "of training data!"))
      } else {
        newdatax <- newdatax[, mgpr$XColName]
      }
    }
    ##############################################
    # Check if normalization is needed for newdata; check model object
    if (is.null(mgpr$trainMeanSd)) { # trainMeanSd includes mean/sd from train.
      normpred <- FALSE
    } else {
      normpred <- TRUE
    }
    
    if (normpred) {
      # Test if training data normalized; stop if not.
      
      # Just logical checks
      if (is.null(mgpr$trainMeanSd)) {
        stop("The training data is not normalized.")
      }
      if (is.null(dim(newdatax))) {
        stop("Error: newdatax is null.")
      }
      
      preds2 <- apply(newdatax, 1, FUN = function(x) (
        (x - mgpr$trainMeanSd[1, ]) / mgpr$trainMeanSd[2, ]))
      
      if (is.null(dim(preds2))) {
        preds2 <- t(as.array(preds2))
      }
      
      newdatax <- NULL
      ##############################################
      Targetdatax <- preds2
    } else {
      if (!is.null(dim(newdatax))) { # Do not transpose vector (i.e. n = 1 plot)
        Targetdatax <- t(newdatax)
      } else {
        Targetdatax <- newdatax
      }
    }
    ##############################################
  }
  # Further argument checks:
  
  # Check CI argument
  if (!is.null(credinter) && !is.numeric(credinter)) {
    stop("Error: Invalid credinter argument.")
  } else {
    if (is.numeric(credinter)) {
      if (credinter <= 0 || credinter >= 1) {
        stop(paste0(
          "Error: Define credible intervals! ",
          "For example set credinter = 0.95 to get 95% credible intervals"
        ))
      }
    }
  }
  
  # Get kernel parameters
  corlen <- mgpr$corlen # correlation length
  ksigma <- mgpr$sigma
  
  # Get noise variance
  errorvar <- mgpr$errorvar
  
  # Check the number of predictors
  # If not vector the upper most and if vector the bottommost
  if (!is.null(dim(Targetdatax))) {
    if (dim(Traindatax)[1] != dim(Targetdatax)[1]) {
      stop("Different number of predictors in train and target data")
    }
  } else {
    if (dim(Traindatax)[1] != 1) {
      stop("Different number of predictors in train and target data")
    }
  }
  
  ##############################################
  n_train <- dim(Traindatax)[2] # number of plots
  
  # If the target is a single plot then:
  if (is.null(dim(Targetdatax))) {
    n_target <- 1 # number of targets
    mn <- length(Targetdatax) # number of predictors
    # Otherwise:
  } else {
    n_target <- dim(Targetdatax)[2] # number of targets
    mn <- dim(Targetdatax)[1] # number of predictors
  }
  
  nx <- dim(t(mgpr$trainy))[1] # number of attributes
  
  ##############################################
  # Objective and gradient for optimization
  trungprobj <- function(x) {
    # negative log-likelihood
    p <- sum((iLx %*% (x - mx))^2)
    return(p)
  }
  trungprgrad <- function(x) {
    # gradient required
    dp <- 2 * t(iLx) %*% (iLx %*% (x - mx))
    return(dp)
  }
  
  ##############################################
  # Matrices used in the prediction.
  
  if (!is.null(dim(Targetdatax))) {
    d <- t(apply(Targetdatax, 2, 
                 FUN = function(x) sqrt(colSums((Traindatax - x)^2))))
  } else { # only one predictor
    d <- t(apply(t(as.array(Targetdatax)), 2, 
                 FUN = function(x) sqrt(colSums((Traindatax - x)^2))))
  }
  Kstar <- mgpr$covfunc(d, ksigma, corlen)
  if (!is.null(dim(Targetdatax))) {
    mu_Target <- apply(Targetdatax, 2, FUN = mgpr$meanf)
  } else { # only one predictor
    mu_Target <- apply(t(as.array(Targetdatax)), 2, FUN = mgpr$meanf)
  }
  
  # univariate case allows simplifications
  if (nx == 1) {
    M3 <- as.vector(mgpr$Cy) * Kstar
    meanpreds_tar <- mu_Target + t(M3 %*% mgpr$predM2)
  } else {
    M3 <- kronecker(mgpr$Cy, Kstar)
    meanpreds_tar <- mu_Target + t(array(M3 %*% 
                                           mgpr$predM2, dim = c(n_target, nx)))
  }
  
  if (!is.null(credinter) | fixneg | covout) {
    ##############################################
    # prediction covariances go here
    covpreds_tar <- array(0, dim = c(nx, nx, n_target))
    # credible intervals go here
    lowerlimgp <- array(0, dim = c(nx, n_target))
    upperlimgp <- array(0, dim = c(nx, n_target))
    
    # sigma value for the credible intervals
    cisig <- qnorm(credinter + (1 - credinter) / 2)
    
    Kstarstar <- mgpr$covfunc(0, ksigma, corlen) * mgpr$Cy + mgpr$E
    for (targetplot in 1:n_target) {
      # indices to target plot rows
      inds <- (seq(targetplot, nx * n_target, n_target))
      if (nx == 1) {
        covpreds_tar[, , targetplot] <- Kstarstar - sum(M3[inds, ] *
                                                          (mgpr$predM1 %*% M3[inds, ]))
      } else {
        covpreds_tar[, , targetplot] <- Kstarstar - M3[inds, ] %*%
          tcrossprod(mgpr$predM1, M3[inds, ])
      }
      
      if (!is.null(credinter)) {
        # Credible intervals
        # special case if response is not multivariate
        if (nx == 1) {
          upperlimgp[, targetplot] <- meanpreds_tar[, targetplot] +
            cisig * sqrt(covpreds_tar[, , targetplot])
          lowerlimgp[, targetplot] <- meanpreds_tar[, targetplot] -
            cisig * sqrt(covpreds_tar[, , targetplot])
        } else {
          upperlimgp[, targetplot] <- meanpreds_tar[, targetplot] +
            cisig * sqrt(diag(covpreds_tar[, , targetplot]))
          lowerlimgp[, targetplot] <- meanpreds_tar[, targetplot] -
            cisig * sqrt(diag(covpreds_tar[, , targetplot]))
        }
        
        ##############################################
        # Upper and lower bounds for fixed CI
        # negative log-likelihood
        # Function trungprobj; read in the beginning of the code.
        if (fixneg && any(lowerlimgp[, targetplot] < 0)) {
          # fix negative CIs
          # check the attributes that need to be iterated
          cond <- lowerlimgp[, targetplot] < 0
          lowerlimgp[cond, targetplot] <- 0
          for (itru in (1:nx)[cond]) {
            # special case if response is not multivariate
            if (nx == 1) {
              wtr <- pnorm(
                0, meanpreds_tar[, targetplot][itru],
                sqrt(covpreds_tar[, , targetplot][1])
              )
              upperlimgp[itru, targetplot] <- qnorm(
                credinter + (1 - credinter) *
                  wtr, meanpreds_tar[, targetplot][itru],
                sqrt(covpreds_tar[, , targetplot][1])
              )
            } else {
              wtr <- pnorm(
                0, meanpreds_tar[, targetplot][itru],
                sqrt(covpreds_tar[, , targetplot][itru, itru])
              )
              upperlimgp[itru, targetplot] <- qnorm(
                credinter + (1 - credinter) *
                  wtr, meanpreds_tar[, targetplot][itru],
                sqrt(covpreds_tar[, , targetplot][itru, itru])
              )
            }
          }
        }
      }
      ##############################################
      # fix negative attribute predictions
      if (fixneg && any(meanpreds_tar[, targetplot] < 0)) {
        btrun <- Inf * array(rep(1, nx), dim = c(1, nx))
        # define pars used in the trungprobj func; Defined outside the function
        fooC <- covpreds_tar[, , targetplot]
        ifooC <- solve(fooC, diag(nx))
        # cholesky of precision matrix
        iLtc <- chol(ifooC)
        # Variables for trungprobj functions
        mx <- meanpreds_tar[, targetplot]
        iLx <- iLtc
        # optimization does not accept initial values x0 < lb
        # replace negs with 0.1
        x00 <- meanpreds_tar[, targetplot]
        x00[x00 < 0] <- 0.1
        # optimize using L-BFGS-B
        meanpreds_tar[, targetplot] <- optim(
          par = x00, fn = trungprobj,
          gr = trungprgrad,
          lower = array(rep(0, nx),
                        dim = c(1, nx)
          ),
          upper = btrun,
          method = "L-BFGS-B"
        )$par
      }
    }
  }
  
  meanpreds_tar <- as.data.frame(t(as.data.frame(meanpreds_tar)))
  
  # Check if column and row names are used.
  if (mgpr$YColName[1] != -9999) {
    colnames(meanpreds_tar) <- mgpr$YColName
  }
  # Rows cannot be added when predicting to new data
  if (mgpr$XRowName[1] != -9999 & (identical(Traindatax, Targetdatax))) {
    rownames(meanpreds_tar) <- mgpr$XRowName
  }
  # For newdatax prediction matrix take the rownames that were saved earlier.
  if (newdatax.rownames[1] != -9999) {
    rownames(meanpreds_tar) <- newdatax.rownames
  }
  if (covout) {
    covpreds_tar <- provideDimnames(covpreds_tar,
                                    base = list(
                                      colnames(meanpreds_tar),
                                      colnames(meanpreds_tar),
                                      paste(rep(
                                        "plot_",
                                        dim(covpreds_tar)[3]
                                      ),
                                      seq(1, dim(covpreds_tar)[3], 1),
                                      sep = ""
                                      )
                                    )
    )
  }
  # data frame out
  if (!is.null(credinter)) {
    CredInter <- data.frame(t(lowerlimgp), t(upperlimgp))
    colnames(CredInter) <- c(
      paste(colnames(meanpreds_tar),
            "_lowlim",
            sep = ""
      ),
      paste(colnames(meanpreds_tar), "_upplim", sep = "")
    )
    out <- list(pred = meanpreds_tar, credinter = CredInter)
    if (covout) {
      out[["covpred"]] <- covpreds_tar
    }
    return(out)
  } else {
    if (covout) {
      out <- list(pred = meanpreds_tar, covpred = covpreds_tar)
    } else {
      out <- as.data.frame(meanpreds_tar)
    }
    return(out)
  }
}