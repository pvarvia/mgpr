mgpr <img src="figs/mgprlogo_nobg.png" width="220" align="right"/>
=================================================
![license](https://img.shields.io/badge/Licence-GPL--3-blue.svg) 

R package for multivariate Gaussian process regression 

To cite the package use `citation()` from within R:

```r
citation("mgpr")
# Varvia, P., Räty, J. and Packalen, P. (2022). mgpr....
```   

# Key functionalities

### Demo data

The *mgprdata* is used for the demonstration of area-based forest inventory with the mgpr library.
The *mgprdata* comprises 825 circular field plots located in the Finnish boreal forests. 
For each field plot, there are 36 plot-level features calculated from the airborne laser scanning (ALS) data.

The *mgprdata* is a part of the forest inventory data used for remote sensing-based forest management inventories in Finland. 
The acquisition of the field data (kaukokartoituskoealat/inventointikoealat) is operated by Finnish forest centre.

The field data belong to [the open data of the Finnish forest centre](https://www.metsakeskus.fi/fi/avoin-metsa-ja-luontotieto/metsatietoaineistot/metsavaratiedot).

The low-density ALS data belong to [the open data of the National Land Survey of Finland](https://www.maanmittauslaitos.fi/en/maps-and-spatial-data/expert-users/product-descriptions/laser-scanning-data-05-p). 

```r
data(mgprdata)
```  

### Fitting an mgpr model

The *mgpr* function is used to fit a Multivariate Gaussian Process Regression model.

```r
gp0 <- mgpr(datay = mgprdata[, 2:5], datax = mgprdata[, 5:29],
kernel = "matern32", kernpar = list(sigma = 1, corlen = 5, errorvar = 0.1))
```   

### Predict method for an mgpr model
The *predict* function is used to predict using an mgpr model. 
The *predict* function supports k-fold cross validation, new predictor variables, limiting predictions to positive scale, and generation of credible intervals.
```r
gp0_pred <- predict(gp0, credinter = 0.95)
```   

### Summarizing an mgpr model

The *summary* function for class "mgpr".

```r
summary(gp0)
```  

# BLAS/LAPACK
The core functionalities of "mgpr" use highly vectorized operations. Thus, we recommend the use of high-performance linear algebra libraries (e.g. Intel oneMKL) linked to R. ...

# Related publications

Varvia, P., Lähivaara, T., Maltamo, M., Packalen, P. and Seppänen, A. (2019). Gaussian Process Regression for Forest Attribute 
Estimation From Airborne Laser Scanning Data, IEEE Transactions on Geoscience and Remote Sensing, vol. 57, no. 6, pp. 3361--3369. https://doi.org/10.1109/TGRS.2018.2883495

Varvia, P., Räty, J., Korhonen, L., and Packalen, P. (2021). Gaussian Process Regression for Airborne Laser Scanning Based Forest Inventory: 
Validation and Parameter Selection. In Proceedings of the SilviLaser Conference 2021, pp. 98--100. https://doi.org/10.34726/wim.1928

Räty, J., Varvia, P., Korhonen, L., Savolainen, P., Maltamo, M. and Packalen, P. (2022). A Comparison of Linear-Mode and Single-Photon Airborne LiDAR in Species-Specific Forest Inventories, 
IEEE Transactions on Geoscience and Remote Sensing, vol. 60, pp. 1--14.  https://doi.org/10.1109/TGRS.2021.3060670

Varvia, P., Räty, J. and Packalen, P. (2022) mgpr: An R ......, TBA 



