test_that("mgpr object is correct (univariate)", {
  #dummy data
  x <- seq(1,10)
  y <- x
  m <- mgpr(datay = y, datax = x,kernel = "rbf",
            kernpar = list(sigma = 1, corlen = 2, errorvar = 0.1))
  
  expect_s3_class(m, "mgpr")
  expect_setequal(m$trainy,y)
  #x values are standardized
  expect_setequal(m$trainx,(x-mean(x))/sd(x))
  #matrix dimensions
  expect_equal(nrow(m$Cy),1)
  expect_equal(ncol(m$Cy),1)
  expect_equal(dim(m$E),NULL) #scalar
  expect_equal(nrow(m$K),length(x))
  expect_equal(ncol(m$K),length(x))
  expect_equal(nrow(m$Ie),length(x))
  expect_equal(ncol(m$Ie),length(x))
  expect_equal(nrow(m$predM1),length(x))
  expect_equal(ncol(m$predM1),length(x))
  expect_equal(nrow(m$predM2),length(x))
  expect_equal(ncol(m$predM2),1)
  #K is positive definite
  expect_no_error(chol(m$K))
  #Cy is positive definite
  expect_no_error(chol(m$Cy))
  #K diagonal elements should be sigma^2
  expect_setequal(diag(m$K),1)
})

test_that("mgpr object is correct (multivariate)", {
  #dummy data
  x <- data.frame(seq(1,10),seq(1,10)^2)
  y <- x
  m <- mgpr(datay = y, datax = x,kernel = "rbf",
            kernpar = list(sigma = 1, corlen = 2, errorvar = 0.1))
  
  expect_s3_class(m, "mgpr")
  expect_setequal(m$trainy[,1],y[,1])
  expect_setequal(m$trainy[,2],y[,2])
  #x values are standardized
  expect_setequal(m$trainx[,1],(x[,1]-mean(x[,1]))/sd(x[,1]))
  expect_setequal(m$trainx[,2],(x[,2]-mean(x[,2]))/sd(x[,2]))
  #matrix dimensions
  expect_equal(nrow(m$Cy),2)
  expect_equal(ncol(m$Cy),2)
  expect_equal(nrow(m$E),2)
  expect_equal(ncol(m$E),2)
  expect_equal(nrow(m$K),nrow(x))
  expect_equal(ncol(m$K),nrow(x))
  expect_equal(nrow(m$Ie),nrow(x))
  expect_equal(ncol(m$Ie),nrow(x))
  expect_equal(nrow(m$predM1),2*nrow(x))
  expect_equal(ncol(m$predM1),2*nrow(x))
  expect_equal(nrow(m$predM2),2*nrow(x))
  expect_equal(ncol(m$predM2),1)
  #K is positive definite
  expect_no_error(chol(m$K))
  #Cy is positive definite
  expect_no_error(chol(m$Cy))
  #K diagonal elements should be sigma^2
  expect_setequal(diag(m$K),1)
})