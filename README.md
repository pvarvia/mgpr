mgpr <img src="figs/mgprlogo_nobg.png" width="220" align="right"/>
=================================================
![license](https://img.shields.io/badge/Licence-GPL--3-blue.svg)

An R package for multivariate Gaussian process regression

To cite the package use `citation()` from within R:

```r
citation("mgpr")
# Varvia, P., Räty, J. and Packalen, P. (2023). mgpr....
```   
# Installation
Install the *mgpr* from GitHub github using the R package *remotes*:
```r
remotes::install_github("pvarvia/mgpr", ref = "main")
```

# Key functionalities

### Fitting an mgpr model

The *mgpr* function is used to fit a multivariate Gaussian process regression model.
By default *mgpr* fits a GP with the Matérn 3/2 kernel/covariance function and training data average mean function. The kernel parameters (length scale, kernel variance and error variance) are optimized using the training data.
```r
gp <- mgpr(datay = mgprdata[, 1:3], datax = mgprdata[, 4:39])
```
The implemented covariance functions include Matérn 1/2, 3/2 (default), and 5/2, and squared exponential/rbf. Mean functions include zero, training data average (default) and a linear trend. These can be specified using the options *kernel* and *meanf*.
```r
gp <- mgpr(datay = mgprdata[, 1:3], datax = mgprdata[, 4:39], kernel = "rbf", meanf = "linear")
```
The kernel parameters can also be set manually using *kernpar*, for example:
```r
gp <- mgpr(datay = mgprdata[, 1:3], datax = mgprdata[, 4:39], 
           kernpar = list(sigma = 1, corlen = 5, errorvar = 0.1))
```
If only some of the kernel parameters are set, the missing parameters will be optimized using training data. For example
```r
gp <- mgpr(datay = mgprdata[, 1:3], datax = mgprdata[, 4:39], kernpar = list(corlen = 5))
```
will optimize kernel variance *sigma* and error variance *errorvar*, while keeping the length scale at 5.

The kernel parameter optimization uses bounded simulated annealing as implemented in the package [optimization](https://cran.r-project.org/web/packages/optimization/index.html). The bounds (*optlower* and *optupper*), starting point (*optstart*), number of data folds (*optkfold*), and the control parameters of the simulated annealing (*optcontrol*) can be modified using *optimpar*, the default values are
```r
gp <- mgpr(datay = mgprdata[, 1:3], datax = mgprdata[, 4:39], 
           optimpar = list(optkfold = 5, optstart = c(1, 10, 0.1),
                           optlower = c(0.3, 3, 0.03), optupper = c(10, 50, 0.5), 
                           optcontrol = list(t0 = 10, nlimit = 50, r = 0.9)))
```
Note that the *optcontrol* list can be used to overwrite the default parameters defined in the [optimization::optim_sa (ver. 1.0-9)](https://cran.r-project.org/web/packages/optimization/optimization.pdf) function.

### Predict method for an mgpr model
The *predict* function is used to predict using a mgpr model. For example, let's split the demo data to separate train and test sets by taking every third row and train a default *mgpr* model using the training data.
```r
dtest <- mgprdata[seq(1, nrow(mgprdata), 3),]
dtrain <- mgprdata[-seq(1, nrow(mgprdata), 3),]
gp <- mgpr(datay = dtrain[, 1:3], datax = dtrain[, 4:39])
```   
We can now predict the response variables on the test data using *predict*:
```r
gp_pred <- predict(gp, newdatax = dtest[, 4:39])
```
This will output a data frame with the predicted values. Usually we also want credible intervals (CI) for the predictions, to get e.g. 95% CIs:
```r
gp_pred <- predict(gp, newdatax = dtest[, 4:39], credinter = 0.95)
```
Setting *covout* will output the prediction covariance matrices:
```r
gp_pred <- predict(gp, newdatax = dtest[, 4:39], covout = TRUE)
```
With either *credinter* or *covout* specified, the output will be a list of data frames, with separate data frames for predictions, credible intervals, and prediction covariances.

The package implements a simple method to bound predictions to be non-negative, which is useful when the response variables represent non-negative attributes, such as volume or height. This is called using the option *fixneg*, for example:
```r
gp_pred <- predict(gp, newdatax = dtest[, 4:39], fixneg = TRUE)
```

The *predict* function can also be used to do k-fold cross-validation using the training data, this is specified by the *kfold* value. For example, to do 10-fold cross-validation:
 ```r
gp_pred <- predict(gp, kfold = 10)
```
Setting *kfold* to the number of training data will do leave-one-out cross-validation.

### Summarizing an mgpr model

The *summary* function for class "mgpr". Prints information on the training data and GP parameters.

```r
summary(gp)
```  

### Demo data

The *mgprdata* is used for the demonstration of area-based forest inventory with the *mgpr* package.
The *mgprdata* comprises 826 circular field plots located in the Finnish boreal forests.
For each field plot, there are 36 plot-level features calculated from the airborne laser scanning (ALS) data using LAStools ([rapidlasso GmbH](http://rapidlasso.com/LAStools)).

The *mgprdata* is a part of the forest inventory data used for remote sensing-based forest management inventories in Finland.
The acquisition of the field data (kaukokartoituskoealat/inventointikoealat) is operated by Finnish forest centre.

The field data belong to [the open data of the Finnish forest centre](https://www.metsakeskus.fi/fi/avoin-metsa-ja-luontotieto/metsatietoaineistot/metsavaratiedot).

The low-density ALS data belong to [the open data of the National Land Survey of Finland](https://www.maanmittauslaitos.fi/en/maps-and-spatial-data/expert-users/product-descriptions/laser-scanning-data-05-p).

```r
data(mgprdata)
```  

# BLAS/LAPACK
The core functions of the mgpr package use vectorized operations. Thus, we recommend the use of high-performance multi-threaded linear algebra libraries (e.g. [Intel oneAPI MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)) linked to R.

# Related publications

Varvia, P., Lähivaara, T., Maltamo, M., Packalen, P. and Seppänen, A. (2019). Gaussian Process Regression for Forest Attribute
Estimation From Airborne Laser Scanning Data, IEEE Transactions on Geoscience and Remote Sensing, vol. 57, no. 6, pp. 3361-3369. https://doi.org/10.1109/TGRS.2018.2883495

Varvia, P., Räty, J., Korhonen, L., and Packalen, P. (2021). Gaussian Process Regression for Airborne Laser Scanning Based Forest Inventory:
Validation and Parameter Selection. In Proceedings of the SilviLaser Conference 2021, pp. 98-100. https://doi.org/10.34726/wim.1928

Räty, J., Varvia, P., Korhonen, L., Savolainen, P., Maltamo, M. and Packalen, P. (2022). A Comparison of Linear-Mode and Single-Photon Airborne LiDAR in Species-Specific Forest Inventories,
IEEE Transactions on Geoscience and Remote Sensing, vol. 60, pp. 1-14.  https://doi.org/10.1109/TGRS.2021.3060670

Varvia, P., Räty, J. and Packalen, P. (2023) mgpr: An R ......, TBA
