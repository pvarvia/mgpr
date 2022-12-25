# The mgpr library for Multivariate Gaussian Process Regression

# Authors: Varvia, P., Räty, J., Packalen, P. 2022.

# Implementation is based on Varvia et al. 2019 in TGRS:
# Varvia, P., Lähivaara, T., Maltamo, M., Packalen, P. and Seppänen, A. (2019). 
# Gaussian Process Regression for Forest Attribute Estimation From 
# Airborne Laser Scanning Data, 
# IEEE Transactions on Geoscience and Remote Sensing, 
# vol. 57, no. 6, pp. 3361-3369. 
# https://doi.org/10.1109/TGRS.2018.2883495

#########################################

#' Fitting an mgpr model
#'
#' \code{mgpr} is used to fit a Multivariate Gaussian Process Regression model.
#' @param datay a data frame or an object coercible by as.data.frame to a
#' data frame containing the response variable or variables (columns) for each
#' observation (rows). It is recommended to use row and column names.
#' @param datax a data frame or an object coercible by as.data.frame to a
#' data frame containing the predictor variables (columns) for each observation
#' (rows). It is recommended to use row and column names.
#' @param kernel a name of the kernel function.
#' Alternatives: "\code{matern12}", "\code{matern32}" (default), "\code{matern52}", and
#' "\code{rbf}".
#' @param kernpar kernel hyperparameters \code{sigma}, correlation length \code{corlen} and noise variance
#' \code{errorvar} given as a list.
#' Kernel hyperparameters not set by the user are estimated using the training set.
#' @param meanf mean function, choices: "\code{zero}", "\code{avg}" (default), "\code{linear}".
#' Mean function specifies a deterministic trend in the Gaussian process and mostly
#' affects extrapolation outside the training set. "\code{zero}" (default)
#' assumes no trend, "\code{avg}" uses the training set mean, "\code{linear}"
#' fits a linear minimum mean square error model to the training set.
#' @param verbose logical. Report extra information on progress.
#' @param optimpar List of control parameters for the hyperparameter estimation:
#' \code{optkfold} defines a kfold value for the cost function (must be in the range
#' 2 and nrow(data)), \code{optstart} defines an initial vector for the optimization,
#' \code{optlower} defines a lower bound vector, \code{optupper} defines a upper bound vector, 
#' and \code{optcontrol} includes the control parameters of the simulated annealing used 
#' as an optimization approach in the hyperparameter estimation.
#' Initial and bound vectors have the order \code{c(sigma, corlen, errorvar)}.
#' For the control parameters, please refer to the documentation of the \code{optimization::optim_sa} function. 
#' Empty list uses the default parameters
#' defined in the \code{optimization::optim_sa} function (ver. 1.0-9). The default setup of
#' \code{mgpr} modifies three of the \code{optimization::optim_sa} function's default parameters: \code{optcontrol = list(t0 = 10, nlimit = 50, r = 0.9)}. 
#'
#' @param ... additional arguments.
#'
#' @return \code{mgpr} returns an object of class "\code{mgpr}". The function
#' \code{\link{summary}} prints a summary of an \code{mgpr} object and the
#' function \code{\link{predict}} is used to make predictions with it.
#'
#' An object of class "\code{mgpr}" is a list containing at least the following
#' components: \tabular{llllllllllllll}{
#'   \code{trainy} \tab \tab response variables of training data. \cr
#'   \code{trainx} \tab \tab predictor variables of training data. \cr
#'   \code{trainMeanSd} \tab \tab mean and standard deviation of training data. \cr
#'   \code{Cy} \tab \tab sample covariance for response variables. \cr
#'   \code{E} \tab \tab error matrix. \cr
#'   \code{Ie} \tab \tab identity matrix for setting E for all observations. \cr
#'   \code{K} \tab \tab kernel matrix. \cr
#'   \code{kernel} \tab \tab a name of the kernel function. \cr
#'   \code{sigma} \tab \tab value of sigma. \cr
#'   \code{corlen} \tab \tab value of correlation length. \cr
#'   \code{errorvar} \tab \tab errorvar noise variance (a.k.a. nugget parameter). \cr
#'   \code{meanf} \tab \tab name of the mean function \cr
#'   \code{XColName} \tab \tab column names of predictor variables. \cr
#'   \code{YColName} \tab \tab column names of response variable(s). \cr
#'   \code{XRowName} \tab \tab row names of predictor variables. \cr
#' }
#' @references Varvia, P., Lahivaara, T., Maltamo, M., Packalen, P. and Seppanen, A. 2019.
#' Gaussian process regression for forest attribute estimation from airborne laser scanning
#' data. IEEE Transactions on Geoscience and Remote Sensing 57(6): 3361-3369.
#'
#' R package article to appear.
#'
#' @examples
#' data(mgprdata)
#' m <- mgpr(
#'   datay = mgprdata[, 1:3], datax = mgprdata[, 4:39],
#'   kernel = "matern32", kernpar = list(sigma = 1, corlen = 5, errorvar = 0.1)
#' )
#' @import Matrix
#' @import optimization
#' @export
mgpr <- function(datay,
                 datax,
                 kernel = "matern32",
                 kernpar = list(sigma = NULL, corlen = NULL, errorvar = NULL),
                 meanf = "avg",
                 optimpar = list(optkfold   = 5,
                                 optstart   = c(1, 10, 0.1),
                                 optlower   = c(0.3, 3, 0.03),
                                 optupper   = c(10, 50, 0.5),
                                 optcontrol = list(t0 = 10, nlimit = 50, 
                                                   r = 0.9)),
                 verbose = FALSE, ...) {
  normpred <- TRUE # Default: always normalize
  # Catch "developer" arguments; normalization can be switched off (for kfold).
  dev_args <- list(...)
  if (length(dev_args) > 1) {
    stop("Error: Unknown argument.")
  } else if (length(dev_args) == 0) { # no developer arguments
    dev_args <- NULL
  } else { # One developer argument
    if (identical(names(dev_args[1]), "normalize")) {
      normpred <- dev_args[[1]]
    } else {
      stop("Error: Unknown argument.")
    }
  }

  if (!is.logical(normpred)) {
    stop("Error: Invalid developer argument.")
  }
  ##############################################
  # Parameters. Check.
  if (missing(datay) || missing(datax)) {
    stop("Training data missing!")
  } else {
    if (any(is.na(datay)) || any(is.na(datax))) {
      stop("Invalid data! [NA values]")
    }
    if (!is.vector(datay) & !is.vector(datax)) { # both df 
      if (dim(datay)[1] != dim(datax)[1]) {
        stop("Invalid data! [Too few training observations]")
      }
      
      if (dim(datay)[1] == 1 | dim(datax)[1] == 1) {
        stop("Invalid data! [Too few training observations]")
      }

    } else if (!is.vector(datay) & is.vector(datax)) { # y is df, x vector
      if (dim(datay)[1] != length(datax)) {
        stop("Invalid data! [Too few training observations]")
      }
      
      if (dim(datay)[1] == 1 | length(datax) == 1) {
        stop("Invalid data! [Too few training observations]")
      }
      if ((dim(datay)[1] <= dim(datay)[2]) | (length(datax) <= dim(datay)[2])) {
        stop("Invalid data! [Too few training observations]")
      }
    } else if (is.vector(datay) & !is.vector(datax)) { # y vector, x df
      if (length(datay) != dim(datax)[1]) {
        stop("Invalid data! [Too few training observations]")
      }
      
      if (length(datay) == 1 | dim(datax)[1] == 1) {
        stop("Invalid data! [Too few training observations]")
      }
      
    } else { # both vectors
      if (length(datay) == 1 | length(datax) == 1) {
        stop("Invalid data! [Too few training observations]")
      }
    }
  }

  if (!is.logical(verbose)) {
    stop("Invalid verbose argument.")
  }
  
  if (identical(kernel, "matern12") || identical(kernel, "matern32") ||
      identical(kernel, "matern52") || identical(kernel, "rbf")) {
    if (verbose) {
      cat(paste0("Using ", kernel, " kernel function"), fill = TRUE)
    }
  } else {
    stop(paste0(
      "Invalid kernel function. Options:",
      "'matern12', 'matern32', 'matern52' or 'rbf'"
    ))
  }
  
  # Kernel parameters
  if (!is.null(kernpar$corlen)) {
    if (is.numeric(kernpar$corlen) &
        length(kernpar$corlen) == 1 &
        kernpar$corlen > 0) {
      corlen <- kernpar$corlen
    } else {
      stop("Invalid kernel parameters: corlen must be a positive scalar.")
    }
  } else {
    corlen <- NULL
  }
  if (!is.null(kernpar$sigma)) {
    if (is.numeric(kernpar$sigma) &
        length(kernpar$sigma) == 1 &
        kernpar$sigma > 0) {
      ksigma <- kernpar$sigma
    } else {
      stop("Invalid kernel parameters: sigma must be a positive scalar.")
    }
  } else {
    ksigma <- NULL
  }
  if (!is.null(kernpar$errorvar)) {
    if (is.numeric(kernpar$errorvar) &
        length(kernpar$errorvar) == 1 &
        kernpar$errorvar > 0) {
      errorvar <- kernpar$errorvar
    } else {
      stop("Invalid kernel parameters: errorvar must be a positive scalar.")
    }
  } else {
    errorvar <- NULL
  }

  # Parameter optimization settings
  optimpar_user <- optimpar
  optimpar <- list(
    optkfold   = 5,
    optstart   = c(1, 10, 0.1),
    optlower   = c(0.3, 3, 0.03),
    optupper   = c(10, 50, 0.5),
    optcontrol = list(t0     = 10, # Our default t0, nlimit and r
                      nlimit = 50,
                      r      = 0.9))
  
  # These are optim_sa default control params (optimization package v. 1.0-9)
  optcontrol_orig = list(vf = NULL, 
                    rf      = 1,     
                    dyn_rf  = TRUE,
                    t0      = 1000,
                    nlimit  = 100,
                    r       = 0.6,
                    k       = 1,
                    t_min   = 0.1,
                    maxgood = 100,
                    stopac  = 30,
                    ac_acc  = NA)
  
  # Check for possible typos in control params
  if (any(!(names(optimpar_user$optcontrol) %in% 
            names(optcontrol_orig)))) {
    stop(paste0("Error: Wrong names detected in the ", 
               "parameter names of optimpar$optcontrol."))
  }
  
  optimpar[names(optimpar_user)] <- optimpar_user
  
  # some checks
  if (!is.list(optimpar$optcontrol)) {
    stop("Error: optimpar$optcontrol must be a list.")
  }
  
  if (length(optimpar$optstart) != 3 | length(optimpar$optlower) != 3 | 
      length(optimpar$optupper) != 3) {
    stop(paste0("Error: optimpar$opstart, optimpar$optlower, and ", 
                "optimpar$optupper must be vectors of 3 values."))
  }
  
  # Ensure that the default parameters are as they were listed above:
  # Replace user-defined
  optcontrol_orig[names(optimpar$optcontrol)] <- optimpar$optcontrol
  # replace to optimpar that is the final set of params
  optimpar$optcontrol <- optcontrol_orig
  
  # Transform data. This allows to set a vector input
  if (is.vector(datay)) {
    datay <- array(datay, dim = c(length(datay), 1))
  } else {
    datay <- datay
  }
  if (is.vector(datax)) {
    datax <- array(datax, dim = c(length(datax), 1))
  } else {
    datax <- datax
  }
  ##############################################
  # Collect predictor names (if data frames are used). The names
  # are used in the selection of correct predictors when predicting
  # to newData
  if (is.data.frame(datax)) {
    col_namex <- colnames(datax)
    row_namex <- rownames(datax)
  } else {
    col_namex <- -9999
    row_namex <- -9999
  }
  if (is.data.frame(datay)) {
    col_namey <- colnames(datay)
  } else {
    col_namey <- -9999
  }

  ##############################################
  # Remove cols var() = 0
  # These checks are not needed if kfold is called
  if (normpred) {
    col_var <- as.numeric(apply(datax, MARGIN = 2, FUN = function(x) var(x)))
    rm_cols <- which(col_var == 0)
    if (is.data.frame(datax)) {
      if (length(rm_cols) > 0) {
        rm_name <- names(datax)[rm_cols]
        datax <- datax[, -rm_cols]
        stop(paste0(
          "The following predictor variables have var(x) == 0: ",
          paste(rm_name, collapse = " "), " at column indices ",
          paste(rm_cols, collapse = " "),  
          ". Please remove them from the training data."
        ))
      }
    } else {
      if (length(rm_cols) > 0) {
        datax <- datax[, -rm_cols]
        stop(paste0(
          "The following predictor variable columns have",
          "var(x) == 0: ", paste(rm_cols, collapse = " "), 
          ". Please remove them from the training data."
        ))
      }
    }
  }
  ##############################################
  # Z score normalization for predictor variables
  if (normpred) {
    # Store train mean and sd (for the normalization of newdata)
    mean_sd_train <- array(0, dim = c(2, dim(datax)[2]))
    preds <- NULL
    for (i in 1:ncol(datax)) {
      mean_sd_train[1, i] <- mean(datax[, i]) # mean
      mean_sd_train[2, i] <- sd(datax[, i]) # sd
      preds <- cbind(preds, (datax[, i] - mean_sd_train[1, i]) /
        mean_sd_train[2, i])
    }
    datax <- NULL
  } else {
    mean_sd_train <- NULL
    preds <- datax
    datax <- NULL
  }
  ##############################################
  # Prepare data for the training
  train_X <- t(preds)
  train_Y <- t(datay)

  # sizes
  n_train <- dim(train_X)[2]
  nx <- dim(train_Y)[1]
  ##############################################
  # Initialize covariance function
  
  covfunc <- init_covf(kernel)

  # Distance matrix D
  D <- apply(train_X, 2, FUN = function(x) sqrt(colSums((train_X - x)^2)))

  # Covariance of response variables
  Cy <- cov(t(train_Y))
  Cy <- array(as.numeric(Cy), dim = c(ncol(Cy), ncol(Cy)))

  Ie <- diag(n_train) # identity matrix for setting E for all obs

  ##############################################
  # Return structure initialization
  return_list <- list(
    "trainy" = t(train_Y),
    "trainx" = t(train_X),
    "trainMeanSd" = mean_sd_train,
    "Cy" = Cy,
    "E" = NULL,
    "Ie" = Ie,
    "K" = NULL,
    "covfunc" = covfunc,
    "sigma" = ksigma,
    "corlen" = corlen,
    "errorvar" = errorvar,
    "meanf" = NULL,
    "XColName" = col_namex,
    "YColName" = col_namey,
    "XRowName" = row_namex,
    "predM1" = NULL,
    "predM2" = NULL
  )
  class(return_list) <- c("mgpr")
  
  # Add mean function
  return_list$meanf <- init_meanf(meanf, return_list, verbose)
  # Compute prior expectancy for training data
  mu_Y <- apply(t(return_list$trainx), 2, FUN = return_list$meanf)
  
  ##############################################
  # Hyperparameter optimization
  # the kernel parameters that are set to NULL will be estimated

  if (is.null(ksigma) || is.null(corlen) || is.null(errorvar)) {
    if (verbose) {
      cat("Optimizing hyperparameters...", fill = TRUE)
	  cat("Number of iterations:", fill = TRUE)
	  iter <- 1
	}

    # Cost function for the hyperparameter optimization
    parcostfun <- function(pars) { # pars: [sigma corlen errorvar]
      # Modify mgpr object
      if (is.null(ksigma)) {
        return_list$sigma <- pars[1]
      }
      if (is.null(corlen)) {
        return_list$corlen <- pars[2]
      }
      if (is.null(errorvar)) {
        return_list$errorvar <- pars[3]
      }
      return_list$K <- covfunc(D, pars[1], pars[2])

      if (dim(Cy)[1] == 1) {
        return_list$E <- pars[3] * diag(Cy)
      } else {
        return_list$E <- pars[3] * diag(diag(Cy))
      }
      
      # precomputed parts of prediction equations
      if (nx == 1) {
        return_list$predM1 <- chol2inv(chol(as.vector(return_list$Cy) * 
                                              return_list$K +
                                              return_list$E * return_list$Ie))
        return_list$predM2 <- return_list$predM1 %*% 
                                       as.vector(return_list$trainy - mu_Y)
      } else {
        return_list$predM1 <- chol2inv(chol(kronecker(return_list$Cy, 
                                                      return_list$K) +
                                              kronecker(return_list$E, 
                                                        return_list$Ie)))
        return_list$predM2 <- return_list$predM1 %*% 
                                       as.vector(return_list$trainy - t(mu_Y))
      }
      
      if (verbose) {
        # update iter number
        cat(iter, "\r") 
		flush.console()
		iter <<- iter + 1
      }

      # SSE cost
	  # The folds change in each iteration on purpose.
	  # This essentially aims to approximate the full k-fold CV
	  # using minibatches of k samples and stochastic optimization.
      yhat <- kfold_mgpr(return_list, kfold = optimpar$optkfold)
      cost <- sum((yhat - return_list$trainy)^2)
      return(cost)
    }

    # Simulated annealing
    bestpars <- optim_sa(
      fun = parcostfun,
      start = optimpar$optstart,
      lower = optimpar$optlower,
      upper = optimpar$optupper,
      control = optimpar$optcontrol
    )
    
    # Modify return_list
    if (is.null(ksigma)) {
      return_list$sigma <- bestpars$par[1]
    }
    if (is.null(corlen)) {
      return_list$corlen <- bestpars$par[2]
    }
    if (is.null(errorvar)) {
      return_list$errorvar <- bestpars$par[3]
    }
    
    return_list$K <- covfunc(D, return_list$sigma, return_list$corlen)
    
    if (dim(Cy)[1] == 1 & dim(Cy)[2] == 1) {
      return_list$E <- return_list$errorvar * diag(Cy)
    } else {
      return_list$E <- return_list$errorvar * diag(diag(Cy))
    }

    if (nx == 1) {
      return_list$predM1 <- chol2inv(chol(as.vector(return_list$Cy) * 
                                            return_list$K +
                                            return_list$E * return_list$Ie))
      return_list$predM2 <- return_list$predM1 %*% 
                                     as.vector(return_list$trainy - mu_Y)
    } else {
      return_list$predM1 <- chol2inv(chol(kronecker(return_list$Cy, 
                                                    return_list$K) +
                                            kronecker(return_list$E, 
                                                      return_list$Ie)))
      return_list$predM2 <- return_list$predM1 %*% 
                                     as.vector(return_list$trainy - t(mu_Y))
    }
    
    if (verbose) {
      cat(
        "\n",  
        "Finished optimization\n",
        "Sigma:",return_list$sigma,"\n",
        "Correlation length:",return_list$corlen,"\n",
        "Error variance:",return_list$errorvar,"\n"
      )
    }
  } else { # if user specified all kernel parameters
    # Kernel matrix K
    return_list$K <- covfunc(D, ksigma, corlen)

    # A case Y is not multivariate.
    if (dim(Cy)[1] == 1 & dim(Cy)[2] == 1) {
      return_list$E <- errorvar * diag(Cy)
    } else {
      return_list$E <- errorvar * diag(diag(Cy))
    }
    
    if (nx == 1) {
      return_list$predM1 <- chol2inv(chol(as.vector(return_list$Cy) * 
                                            return_list$K +
                                            return_list$E * return_list$Ie))
      return_list$predM2 <- return_list$predM1 %*% 
                                     as.vector(return_list$trainy - mu_Y)
    } else {
      return_list$predM1 <- chol2inv(chol(kronecker(return_list$Cy, 
                                                    return_list$K) +
                                            kronecker(return_list$E, 
                                                      return_list$Ie)))
      return_list$predM2 <- return_list$predM1 %*% 
                                    as.vector(return_list$trainy - t(mu_Y))
    }
  }
  
  return(return_list)
}

#' Summarizing an mgpr model
#'
#' \code{summary} method for class "\code{mgpr}".
#' @param mgpr a \code{mgpr} model object, usually, a result of a call to \code{\link{mgpr}}.
#' @return A summary of the mgpr model will be printed to the R console.
#' @examples
#' data(mgprdata)
#' m <- mgpr(
#'   datay = mgprdata[, 1:3], datax = mgprdata[, 4:40],
#'   kernel = "matern32", kernpar = list(sigma = 1, corlen = 5, errorvar = 0.1)
#' )
#' summary(m)
#' @export
summary <- function(mgpr) {
  UseMethod("summary")
}

#' @method summary mgpr
#' @export
summary.mgpr <- function(mgpr) {
  cat(
    "Multivariate Gaussian Process Regression model:\n",
    dim(mgpr$trainx)[1], "training observations\n",
    dim(mgpr$trainx)[2],
    " predictor variables\n", dim(mgpr$trainy)[2], " response variables",
    "\n\n", "Kernel parameters:",
    "\n", "Kernel function:", environment(mgpr[["covfunc"]])[["covftype"]],
    "\n", "Sigma: ", mgpr$sigma,
    "\n", "Correlation length: ", mgpr$corlen,
    "\n", "Error variance:", mgpr$errorvar,
    "\n", "Mean function: ", environment(mgpr[["meanf"]])[["meanftype"]], "\n"
  )
}

#' Predict function for an mgpr model
#'
#' Predict using a Multivariate Gaussian process regression (mgpr) model. The
#' function supports k-fold cross validation, new predictor variables, limiting
#' predictions to positive scale, and generation of credible intervals.
#' @param mgpr a \code{mgpr} model object, usually, a result of a call to
#' \code{\link{mgpr}}.
#' @param newdatax a data frame or an object coercible by as.data.frame to a
#' data frame containing new predictor variables. Column names must match with
#' the ones in the \code{mgpr} object. If NULL \code{datax} used in model
#' training is used.
#' @param credinter threshold of the credible intervals, e.g. \code{credinter=0.95}.
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
#'   datay = mgprdata[, 1:3], datax = mgprdata[, 4:40],
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
