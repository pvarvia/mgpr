#####################################################################
#
#  Test script of MGPR.
#
#    March 3, 2021
#    Petteri Packalen
#    Janne R?ty
#
#####################################################################
rm( list = ls() )

# setMKLthreads(2)

require(devtools) # Install the newest if not alreary installed
install_github("JJRaty/MGPR", auth_token = "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", force = TRUE)

require(mgpr)
#
# Test case 1: 100 obs, 3 ys and 25 xs.
#
data(mgprdata)

test_x <-  mgprdata[1:50, 5:29]
test_y <-  mgprdata[1:50, 2:4]

test_mgpr <- function(dy, dx, manustep = TRUE) {
manustep <- TRUE

# Performance assessment function. Inputs are data.frames
# having corresponding columns with matching column names.
# Returns absolute and relative rmse and bias inside of a
# data.frame.
#
validate <- function(pred, obs) {
  if (is.numeric(pred) & is.numeric(obs)){
     pred = as.data.frame( pred )
     obs = as.data.frame( obs )
     colnames(pred) = colnames(obs) = "unknown"
  }

  if ( !is.data.frame(pred) | !is.data.frame(obs)) {
    print("ERROR: pred and obs must be of type 'data.frame'.")
    return( NULL )
  }
  if ( mean(colnames(pred)==colnames(obs)) != 1 ) {
    print("ERROR: colnames in pred and obs in validate do not match.")
    return( NULL )
  }
  m = matrix( 0, nrow=4, ncol=ncol(pred) )

  for( i in 1:ncol(pred)) {
    m[1,i] = sqrt( sum( ( pred[,i] - obs[,i] )^2 ) / nrow(obs) )
    m[2,i] = 100 * m[1,i] / mean( obs[,i] )
    m[3,i] = mean( obs[,i] - pred[,i] )
    m[4,i] = 100 * m[3,i] / mean( obs[,i] )
  }

  colnames(m) = colnames( pred )
  rownames(m) = c( "rmse", "rmse%", "bias", "bias%" )
  m = round( m, 2 )

  return( m )
}

# 
cat("Checking if any var(x) == 0: ", fill = TRUE)
tryCatch({m_t1 <- mgpr(dy, dx)},
         error = function(e) print(e))
omit <- apply(dx, 2, var)
if (any(omit == 0)) {
  cat("Removing x variables with var(x) == 0", fill = TRUE)
  dx <- dx[, -which(omit == 0)]
}

# train & predict
cat("Testing Kernel arguments: ", fill = TRUE)
cat("\nTry with a missing kernel argument: ", fill = TRUE)
tryCatch({m_t1 <- mgpr(dy, dx)},
         error = function(e) print(e))
cat("\nTest with an invalid kernel name: ", fill = TRUE)
tryCatch({m_t1 <- mgpr(datay = dy, datax = dx, kernel = "mate")},
         error = function(e) print(e)) # error
tryCatch({m_t1 <- mgpr(datay = dy, datax = dx, kernel = 1)},
         error = function(e) print(e)) # error

cat("\nTesting kernel parameter arguments: ", fill = TRUE)
cat("No hp tuning, with different kernels (no errors)", fill = TRUE)
tryCatch({m_t1_1 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
             kernpar = list(sigma = 1, corlen = 10, errorvar = .1))}) # Works
tryCatch({m_t1_2 <- mgpr(datay = dy, datax = dx, kernel = "matern12",
                       kernpar = list(sigma = 1, corlen = 10, errorvar = .1))}) # Works
tryCatch({m_t1_3 <- mgpr(datay = dy, datax = dx, kernel = "matern52",
                       kernpar = list(sigma = 1, corlen = 10, errorvar = .1))}) # Works
tryCatch({m_t1_4 <- mgpr(datay = dy, datax = dx, kernel = "rbf",
                       kernpar = list(sigma = 1, corlen = 10, errorvar = .1))}) # Works
tryCatch({summary(m_t1_1)}) 
cat("\n", fill =TRUE)
tryCatch({summary(m_t1_2)})
cat("\n", fill =TRUE)
tryCatch({summary(m_t1_3)})
cat("\n", fill =TRUE)
tryCatch({summary(m_t1_4)})
cat("\n", fill =TRUE)

if (manustep) {
  readline(prompt = "\nDoes it look good? Press [enter] to continue")
}


# Define kern hpars manually
cat("\nTesting kernparn and meanf arguments", 
    fill = TRUE)
cat("Empty list (hp optimization)", fill = TRUE)
cat("\n", fill =TRUE)
tryCatch({m_t2_1 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
             kernpar = list())}) 
tryCatch({summary(m_t2_1)})

cat("\nUser-defined corlen = 1 (optimize rest): ", fill = TRUE)
cat("\n", fill =TRUE)
tryCatch({m_t2_2 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
             kernpar = list(corlen = 1))}) 
tryCatch({summary(m_t2_2)}) 

cat("\nUser-defined corlen = 5 and sigma = 1 (optimize rest): ", fill = TRUE)
cat("\n", fill =TRUE)
tryCatch({m_t2_3 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
             kernpar = list(corlen = 5, sigma = 1))}) # works
 tryCatch({summary(m_t2_3)})

cat("\nInvalid meanf: ", fill = TRUE)
cat("\n", fill =TRUE)
tryCatch({m_t2_4 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
                         kernpar = list(corlen = 5, sigma = 1), meanf = 1)})
tryCatch({m_t2_5 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
                         kernpar = list(corlen = 5, sigma = 1), meanf = "xac")}) 

cat("\nValid mean functions: ", fill = TRUE)
cat("\n", fill =TRUE)
tryCatch({m_t2_6 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
                         kernpar = list(corlen = 5, sigma = 1), meanf = "zero")}) # works
tryCatch({summary(m_t2_6)})
cat("\n", fill =TRUE)
tryCatch({m_t2_7 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
                         kernpar = list(corlen = 5, sigma = 1), meanf = "linear")}) # works
tryCatch({summary(m_t2_7)})

if (manustep) {
  readline(prompt = "\nDoes it look good? Press [enter] to continue")
}

cat("\nTesting prediction function:", fill = TRUE)

tryCatch({pred_m_t2_1_tm <- predict(m_t1_1, verbose = F, 
                                    covout = F)},
         error = function(e) print(e)) 

tryCatch({pred_m_t2_1_kfold <- predict( m_t2_1, verbose = F, kfold = 10,
                                covout = F)},
         error = function(e) print(e)) 

tryCatch({pred_m_t2_1_loocv <- predict( m_t2_1, verbose = F, kfold = nrow(dx),
                                        covout = F)},
         error = function(e) print(e)) 
cat("\nTRAIN mode:", fill = TRUE)
print(names(pred_m_t2_1_tm))
#print(head(pred_m_t2_1_tm$pred));print(head(pred_m_t2_1_tm$credinter))
print(validate(pred_m_t2_1_tm, dy))

cat("\n10-fold mode:", fill = TRUE)
print(names(pred_m_t2_1_kfold))
#print(head(pred_m_t2_1_kfold$pred));print(head(pred_m_t2_1_kfold$credinter))
print(validate(pred_m_t2_1_kfold, dy))

cat("\nLOOCV-fold mode:", fill = TRUE)
print(names(pred_m_t2_1_loocv))
#print(head(pred_m_t2_1_kfold$pred));print(head(pred_m_t2_1_kfold$credinter))
print(validate(pred_m_t2_1_loocv, dy))

if (manustep) {
  readline(prompt = "\nDoes it look good? Press [enter] to continue")
}

# Test fixneg argument
cat("\nPost-process negative predictions (should not give errors): ", 
    fill = TRUE)
dy_neg <- dy
dy_neg[, 1:2] <- dy_neg[, 1:2] - 10
m_t1_1_neg <- mgpr(datay = dy_neg, datax = dx, kernel = "matern32",
               kernpar = list(sigma = 1, corlen = 10, errorvar = .1))
tryCatch({print(min(predict(m_t1_1_neg, kfold = 0, verbose = FALSE, fixneg = FALSE)))},
         error = function(e) print(e)) # NO CV
tryCatch({print(min(predict(m_t1_1_neg, kfold = 0, verbose = FALSE, fixneg = TRUE)))},
         error = function(e) print(e))  # NO CV
tryCatch({print(min(predict(m_t1_1_neg, kfold = 10, verbose = FALSE, fixneg = TRUE)))},
         error = function(e) print(e)) # KFOLD
tryCatch({print(min(predict(m_t1_1_neg, kfold = 10, verbose = FALSE, fixneg = FALSE)))},
         error = function(e) print(e)) #KFOLD


if (manustep) {
  readline(prompt = "\nDoes it look good? Press [enter] to continue")
}

# include also credible intervals
cat("\nCheck credible intervals (should not give errors): ", 
    fill = TRUE)
tryCatch({pred_lcv_ci <- predict(m_t2_1, verbose = F, kfold = 10, credinter = 0.95)},
         error = function(e) print(e))
tryCatch({pred_loocv_ci <- predict(m_t2_1, verbose = F, kfold = nrow(dx), credinter = 0.95)},
         error = function(e) print(e))
cat("Compare results from Kfold and LOOCV: ", fill = TRUE)
print(validate(pred_lcv_ci$pred, dy))
print(validate(pred_loocv_ci$pred, dy))
print(head(pred_lcv_ci$credinter)) # lets see credible intervals
print(head(pred_loocv_ci$credinter)) # lets see credible intervals
cat("\nCheck negative credible intervals (fixneg = FALSE versus fixneg = TRUE): ")
tryCatch(pred_lcv_ci_neg <- predict(m_t1_1_neg, verbose = F, kfold = 10, credinter = 0.95))
tryCatch(pred_lcv_ci_negfix <- predict(m_t1_1_neg, verbose = F, kfold = 10, credinter = 0.95, fixneg = TRUE))
print(head(pred_lcv_ci_neg$credinter)) # lets see credible intervals
print(head(pred_lcv_ci_negfix$credinter)) # lets see credible intervals


if (manustep) {
  readline(prompt = "\nDoes it look good? Press [enter] to continue")
}

# predictions with and without credible intervals should be equal
cat("\nCheck that predictions are similar with and without CI mode:")
if (mean(pred_m_t2_1_loocv == pred_loocv_ci) == 1) {
  cat("Predictions OK.", fill = TRUE)
} else {
  stop("Error. Predictions are dissimilar.")
}

# inputs as vectors; unimode
cat("\nTesting vector input for response (univariate):", fill = TRUE)
tryCatch(m_t1_s <- mgpr( as.numeric(dy[, 1]), dx ))
tryCatch(pred_s1 <- predict( m_t1_s, verbose = F))
print(summary(m_t1_s))
print(validate( pred_s1[,1], dy[,1] ))

if (manustep) {
  readline(prompt = "Does it look good? Press [enter] to continue")
}

# one response and one predictor: train & predict
cat("\nTesting vector input for response and predictor variable" , fill = TRUE)
tryCatch( m_t1_s_s <- mgpr( as.numeric(dy[,1]), dx[,1]) )
tryCatch( pred_s_s1 <- predict( m_t1_s_s, verbose=F ))
tryCatch( pred_s_s2 <- predict( m_t1_s_s, verbose=Fc, kfold = 10))
print(summary( m_t1_s_s1 ))
print(validate( pred_s_s1[,1], dy[,1] ))
print(validate( pred_s_s2[,1], dy[,1] ))


cat("\nTesting if NAs in the response data: " , fill = TRUE)
dy_na <- dy
dy_na[sample(1, nrow(dy), 3), 2] <- NA
tryCatch({m_na_t <- mgpr(datay = dy_na, datax = dx, kernel = "matern32",
                         kernpar = list())},
         error = function(e) print(e)) 

cat("Testing if NAs in the predictor variable data: " , fill = TRUE)
dx_na <- dx
dx_na[sample(1, nrow(dy), 3), 2] <- NA
tryCatch({m_dxna_t <- mgpr(datay = dy, datax = dx_na, kernel = "matern32",
                         kernpar = list())},
         error = function(e) print(e))
}

# Run with the mgprdata
test_mgpr(dy = test_y, dx = test_x, manustep = TRUE)
