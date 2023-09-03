#' Summarizing an mgpr model
#'
#' \code{summary} method for class "\code{mgpr}".
#' @param mgpr a \code{mgpr} model object, usually, 
#' a result of a call to \code{\link{mgpr}}.
#' @return A summary of the mgpr model will be printed to the R console.
#' @examples
#' data(mgprdata)
#' m <- mgpr(
#'   datay = mgprdata[, 1:3], datax = mgprdata[, 4:39],
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