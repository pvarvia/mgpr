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
For each field plot, there are ? plot-level features calculated from the airborne laser scanning (ALS) data.

The *mgprdata* is a part of the forest inventory data used for remote sensing-based forest management inventories in Finland. 
The acquisition of the field data (kaukokartoituskoealat) is operated by Finnish forest centre.

The field data are openly available at: https://www.metsakeskus.fi/fi/avoin-metsa-ja-luontotieto/metsatietoaineistot/metsavaratiedot

The low-density ALS data belongs to the open data of the National Land Survey of Finland: https://www.maanmittauslaitos.fi/en/maps-and-spatial-data/expert-users/product-descriptions/laser-scanning-data-05-p 

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

