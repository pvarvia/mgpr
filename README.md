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

### Fitting an MGPR model

*mgpr* is used to fit a Multivariate Gaussian Process Regression model.

```r
data(mgprdata)
gp0 <- mgpr(datay = mgprdata[, 2:5], datax = mgprdata[, 5:29],
kernel = "matern32", kernpar = list(sigma = 1, corlen = 5, errorvar = 0.1))
```   

### Predict method for an mgpr model
Predict using an mgpr model. 
The function supports k-fold cross validation, new predictor variables, limiting predictions to positive scale, and generation of credible intervals.
```r
gp0_pred <- predict(gp0, credinter = 0.95)
```   

### Summarizing an mgpr model

*summary* method for class "mgpr".

```r
summary(gp0)
```  

### Demo data
```r
data(mgprdata)
```  