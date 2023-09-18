test_that("predict output is ok (univariate)", {
  x <- seq(1, 10)
  y <- x
  xnew <- as.data.frame(seq(0.5, 10.5))
  m <- mgpr(datay = y,
            datax = x,
            kernpar = list(sigma = 1, corlen = 2, errorvar = 0.1))
  p <- predict(m,newdatax = xnew, credinter = 0.5)
  expect_equal(nrow(p$pred), 11)
  expect_equal(nrow(p$credinter), 11)
  expect_equal(ncol(p$credinter), 2)
  #upper CI limits should be larger than lower
  expect_true(all(p$credinter$V1_upplim > p$credinter$V1_lowlim))
})

test_that("predict output is ok (multivariate)", {
  #dummy data
  x <- data.frame(x1 = seq(1, 10), x2 = seq(1, 10)^2)
  y <- x
  xnew <- data.frame(x1 = seq(0.5, 10.5),x2 = seq(0.5, 10.5)^2)
  m <- mgpr(datay = y,
            datax = x,
            kernpar = list(sigma = 1, corlen = 2, errorvar = 0.1))
  p <- predict(m,newdatax = xnew, credinter = 0.5)
  expect_equal(nrow(p$pred), 11)
  expect_equal(ncol(p$pred), 2)
  expect_equal(nrow(p$credinter), 11)
  expect_equal(ncol(p$credinter), 4)
  #upper CI limits should be larger than lower
  expect_true(all(p$credinter$x1_upplim > p$credinter$x1_lowlim))
  expect_true(all(p$credinter$x2_upplim > p$credinter$x2_lowlim))
})

test_that("predicted values are ok (simple test)", {
  x <- c(-1, 1)
  y <- x
  m <- mgpr(datay = y,
            datax = x,
            kernpar = list(sigma = 1, corlen = 2, errorvar = 0.1))
  p <- predict(m,newdatax = as.data.frame(0), credinter = 0.5)
  #should be zero by symmetry for all kernel parametrizations
  expect_equal(as.numeric(p$pred), 0)
})

test_that("k-fold output is ok", {
  #dummy data
  x <- data.frame(x1 = seq(1, 10), x2 = seq(1, 10)^2)
  y <- x
  m <- mgpr(datay = y,
            datax = x,
            kernpar = list(sigma = 1, corlen = 2, errorvar = 0.1))
  p <- predict(m, credinter = 0.5, kfold = 2)
  expect_equal(nrow(p$pred), 10)
  expect_equal(ncol(p$pred), 2)
  expect_equal(nrow(p$credinter), 10)
  expect_equal(ncol(p$credinter), 4)
})

test_that("fixneg works", {
  x <- seq(-5, 5)
  y <- x
  xnew <- seq(-5.5,5.5)
  m <- mgpr(datay = y,
            datax = x,
            kernpar = list(sigma = 1, corlen = 2, errorvar = 0.1))
  p <- predict(m,newdatax = as.data.frame(xnew), credinter = 0.5, fixneg = T)
  expect_true(all(p$pred >= 0))
  expect_true(all(p$credinter$V1_lowlim >= 0))
})

test_that("trained parameters work", {
  x <- seq(1, 10)
  y <- x
  m <- mgpr(datay = y, datax = as.data.frame(x))
  expect_no_error(p <- predict(m))
  expect_equal(mean(p$V1-y), 0, tolerance = 1e-3)
})