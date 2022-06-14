#' Dataset for the demonstration of area-based forest inventory with mgpr 
#' 
#' @description 
#' The mgprdata comprises 825 circular field plots (r = 9 m) located in the Finnish boreal forests. 
#' For each field plot, there are 36 plot-level features calculated from the airborne laser scanning (ALS) data. 
#' @details 
#' Observed variables: \tabular{lll}{
#'   \code{Volume, v} \tab \tab Observed timber volume (m3/ha). \cr
#'   \code{Mean height, h} \tab \tab Mean of tree heights (m). \cr
#'   \code{Mean diameter, d} \tab \tab Mean of tree diameters at breast height (cm). \cr
#' }
#' 
#' Features: \tabular{llllll}{
#'   \code{hstd} \tab \tab Standard deviation of ALS-based height measurements. \cr
#'   \code{hmean} \tab \tab Mean of ALS-based height measurements. \cr
#'   \code{hskew} \tab \tab Skewness of the distribution of the ALS-based height measurements. \cr
#'   \code{h10, h30, ..., h90} \tab \tab Percentiles calculated based on the distribution of the ALS-based height measurements.  \cr
#'   \code{dX} \tab \tab Density features calculated based on the distribution of the ALS-based height measurements with a fixed X height threshold.
#'   For example, f_d2 means the proportion of first echoes below 2 m to all first echoes.  \cr
#' }
#' 
#' The features were calculated separately for the last and first echoes (f and l).
#' 
#' @encoding UTF-8
#' @references Varvia, P., Räty, J., Packalen, P. 2022.
#' Paper to be published in 2022.
#' 
#' The mgprdata is a part of the forest inventory data used for remote sensing-based forest inventories in Finland. 
#' The acquisition of the field data (kaukokartoituskoealat) is operated by Finnish forest centre. \cr\cr
#' The field data are openly available at: https://www.metsakeskus.fi/fi/avoin-metsa-ja-luontotieto/metsatietoaineistot/metsavaratiedot
#' \cr\cr
#' The low-density ALS data provided by the National Land Survey of Finland are openly available at: https://www.maanmittauslaitos.fi/en/maps-and-spatial-data/expert-users/product-descriptions/laser-scanning-data-05-p
#'
#' @examples
#' data(mgprdata)
#' 
#' # Univariate model using the matern 32 kernel
#' gp1 <- mgpr(datay = mgprdata$h, datax = mgprdata$f_h90, kernel = "matern32",
#' meanf="avg", kernpar = list(sigma = 1, corlen = 5, errorvar = 0.1))
#' 
#' #' # Multivariate model using the matern 32 kernel
#' gp2 <- mgpr(datay = mgprdata[, c("h", "d", "v")], 
#' datax = mgprdata[, 7:42], kernel = "rbf", meanf = "avg")
#' )
#' @name mgprdata
#' @docType data
#' @author 
#' @usage data(mgprdata)
NULL