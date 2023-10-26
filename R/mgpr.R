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
#' Alternatives: "\code{matern12}", "\code{matern32}" (default), 
#' "\code{matern52}", and
#' "\code{rbf}".
#' @param kernpar kernel hyperparameters \code{sigma}, correlation 
#' length \code{corlen} and noise variance
#' \code{errorvar} given as a list.
#' Kernel hyperparameters not set by the user are 
#' estimated using the training set.
#' @param meanf mean function, choices: "\code{zero}", "\code{avg}" (default),
#'  "\code{linear}".
#' Mean function specifies a deterministic trend in the Gaussian process and mostly
#' affects extrapolation outside the training set. "\code{zero}" (default)
#' assumes no trend, "\code{avg}" uses the training set mean, "\code{linear}"
#' fits a linear minimum mean square error model to the training set.
#' @param verbose logical. Report extra information on progress.
#' @param optimpar List of control parameters for the hyperparameter estimation:
#' \code{optkfold} defines a kfold value for the cost function 
#' (must be in the range 2 and nrow(data)), \code{optstart} defines an initial
#'  vector for the optimization,
#' \code{optlower} defines a lower bound vector, \code{optupper} defines 
#' a upper bound vector, 
#' and \code{optcontrol} includes the control parameters of the simulated 
#' annealing used 
#' as an optimization approach in the hyperparameter estimation.
#' Initial and bound vectors have the order \code{c(sigma, corlen, errorvar)}.
#' For the control parameters, please refer to the documentation of the 
#' \code{optimization::optim_sa} function. 
#' Empty list uses the default parameters
#' defined in the \code{optimization::optim_sa} function (ver. 1.0-9). 
#' The default setup of
#' \code{mgpr} modifies three of the \code{optimization::optim_sa} function's 
#' default parameters: \code{optcontrol = list(t0 = 10, nlimit = 50, r = 0.9)}. 
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
#'   \code{trainMeanSd} \tab \tab mean and standard deviation of 
#'   training data. \cr
#'   \code{Cy} \tab \tab sample covariance for response variables. \cr
#'   \code{E} \tab \tab error matrix. \cr
#'   \code{Ie} \tab \tab identity matrix for setting E for all observations. \cr
#'   \code{K} \tab \tab kernel matrix. \cr
#'   \code{kernel} \tab \tab a name of the kernel function. \cr
#'   \code{sigma} \tab \tab value of sigma. \cr
#'   \code{corlen} \tab \tab value of correlation length. \cr
#'   \code{errorvar} \tab \tab errorvar noise variance 
#'   (a.k.a. nugget parameter). \cr
#'   \code{meanf} \tab \tab name of the mean function \cr
#'   \code{XColName} \tab \tab column names of predictor variables. \cr
#'   \code{YColName} \tab \tab column names of response variable(s). \cr
#'   \code{XRowName} \tab \tab row names of predictor variables. \cr
#' }
#' @references Varvia, P., Lähivaara, T., Maltamo, M., Packalen, P. 
#' and Seppänen, A. (2019). Gaussian process regression for forest attribute 
#' estimation from airborne laser scanning data. 
#' IEEE Transactions on Geoscience and Remote Sensing 57(6): 3361-3369.
#' https://doi.org/10.1109/TGRS.2018.2883495
#'
#' Varvia, P., Räty, J. and Packalen, P. (2023). mgpr: An R package for
#' multivariate Gaussian process regression, SoftwareX 24: 101563.
#' https://doi.org/10.1016/j.softx.2023.101563
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
