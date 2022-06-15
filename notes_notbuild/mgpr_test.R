#####################################################################
# The mgpr library, Varvia, Räty and Packalen 2020
#
#  The test function of the mgpr library
#
#####################################################################
rm(list = ls())

# Install the mgpr library
# require(devtools) # Install the newest if not alreary installed
# install_github("JJRaty/MGPR", 
# ***REMOVED***, force = TRUE)
#
require(mgpr)
# Load data
data(mgprdata)

test_x <-  mgprdata[1:50, 5:40]
test_y <-  mgprdata[1:50, 2:4]

# Run example
# Inputs must be data frames with several observations, several Xs and Ys
# test_mgpr(dy = test_y, dx = test_x, manustep = TRUE)

# The test function
test_mgpr <- function(dy, dx, manustep = TRUE) {
  
  # Performance assessment function.
  # Written by P. Packalen
  validate <- function(pred, obs) {
    if (is.numeric(pred) & is.numeric(obs)){
       pred <- as.data.frame( pred )
       obs <- as.data.frame( obs )
       colnames(pred) <- colnames(obs) <- "unknown"
    }
  
    if (!is.data.frame(pred) | !is.data.frame(obs)) {
      print("ERROR: pred and obs must be of type 'data.frame'.")
      return( NULL )
    }
    if (mean(colnames(pred) == colnames(obs)) != 1) {
      print("ERROR: colnames in pred and obs in validate do not match.")
      return( NULL )
    }
    m <- matrix(0, nrow = 4, ncol = ncol(pred))
  
    for (i in 1:ncol(pred)) {
      m[1,i] <- sqrt( sum((pred[, i] - obs[, i])^2) / nrow(obs))
      m[2,i] <- 100 * m[1, i] / mean(obs[, i])
      m[3,i] <- mean(obs[, i] - pred[, i])
      m[4,i] <- 100 * m[3, i] / mean(obs[, i])
    }
  
    colnames(m) <- colnames(pred)
    rownames(m) <- c("rmse", "rmse%", "bias", "bias%")
    m <- round(m, 2)
  
    return(m)
  }
  
  # 
  cat("\n#################################", fill = TRUE)
  cat("Checking if any var(x) == 0: (one was generated)", fill = TRUE)
  dx$generated_var0 <- 1
  tryCatch({m_t1 <- mgpr(dy, dx)},
           error = function(e) print(e))
  omit <- apply(dx, 2, var)
  
  if (any(omit == 0)) {
    cat(paste0("The test program removed the x variables", 
               " with var(x) == 0. Let's continue."), fill = TRUE)
    dx <- dx[, -which(omit == 0)]
  }
  cat("\n#################################", fill = TRUE)
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
  cat(paste0("No parameter optimization,", 
             "with different kernels (should not return errors)", fill = TRUE))
  tryCatch({m_t1_1 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
               kernpar = list(sigma = 1, corlen = 10, errorvar = .1))}) # Works
  tryCatch({m_t1_2 <- mgpr(datay = dy, datax = dx, kernel = "matern12",
                         kernpar = list(sigma = 1, 
                                        corlen = 10, errorvar = .1))}) # Works
  tryCatch({m_t1_3 <- mgpr(datay = dy, datax = dx, kernel = "matern52",
                         kernpar = list(sigma = 1, 
                                        corlen = 10, errorvar = .1))}) # Works
  tryCatch({m_t1_4 <- mgpr(datay = dy, datax = dx, kernel = "rbf",
                         kernpar = list(sigma = 1, 
                                        corlen = 10, errorvar = .1))}) # Works
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
  
  cat("\n#################################", fill = TRUE)
  # Define kern hpars manually
  cat("\nTesting the kernparn and meanf arguments", fill = TRUE)
  cat("Empty list given and the mgpr library uses the parameter optimization)",
      fill = TRUE)
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
                           kernpar = list(corlen = 5, sigma = 1), 
                           meanf = 1)},
           error = function(e) print(e))
  tryCatch({m_t2_5 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
                           kernpar = list(corlen = 5, sigma = 1), 
                           meanf = "xac")},
           error = function(e) print(e)) 
  
  cat("\nValid mean functions: ", fill = TRUE)
  cat("\n", fill =TRUE)
  tryCatch({m_t2_6 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
                           kernpar = list(corlen = 5, sigma = 1), 
                           meanf = "zero")},
           error = function(e) print(e)) # works
  tryCatch({summary(m_t2_6)})
  cat("\n", fill =TRUE)
  tryCatch({m_t2_7 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
                           kernpar = list(corlen = 5, sigma = 1), 
                           meanf = "linear")},
           error = function(e) print(e)) # works
  tryCatch({summary(m_t2_7)},
           error = function(e) print(e))
  
  if (manustep) {
    readline(prompt = "\nDoes it look good? Press [enter] to continue")
  }
  
  cat("\n#################################", fill = TRUE)
  
  cat(paste0("\nTesting the predict function. " ,
             "Check manually that error rates are reasonable:"),
      fill = TRUE)
  
  tryCatch({pred_m_t1_1_tm <- predict(m_t1_1, verbose = F, covout = F)},
           error = function(e) print(e)) 
  
  tryCatch({pred_m_t1_1_kfold <- predict(m_t1_1, verbose = F, kfold = 10,
                                  covout = F)},
           error = function(e) print(e)) 
  
  tryCatch({pred_m_t1_1_loocv <- predict(m_t1_1, verbose = F, kfold = nrow(dx),
                                          covout = F)},
           error = function(e) print(e)) 
  
  tryCatch({pred_m_t1_1_new <- predict(m_t1_1, verbose = F, newdatax = dx,
                                          covout = F)},
           error = function(e) print(e)) 
  cat("\nTRAIN mode:", fill = TRUE)
  print(names(pred_m_t1_1_tm))
  #print(head(pred_m_t2_1_tm$pred));print(head(pred_m_t2_1_tm$credinter))
  print(validate(pred_m_t1_1_tm, dy))
  
  cat("\n10-fold mode:", fill = TRUE)
  print(names(pred_m_t1_1_kfold))
  #print(head(pred_m_t2_1_kfold$pred));print(head(pred_m_t2_1_kfold$credinter))
  print(validate(pred_m_t1_1_kfold, dy))
  
  cat("\nLOOCV-fold mode:", fill = TRUE)
  print(names(pred_m_t1_1_loocv))
  #print(head(pred_m_t2_1_kfold$pred));print(head(pred_m_t2_1_kfold$credinter))
  print(validate(pred_m_t1_1_loocv, dy))
  
  cat("\nNew data mode:", fill = TRUE)
  print(names(pred_m_t1_1_new))
  #print(head(pred_m_t2_1_kfold$pred));print(head(pred_m_t2_1_kfold$credinter))
  print(validate(pred_m_t1_1_new, dy))
  
  if (manustep) {
    readline(prompt = "\nDoes it look good? Press [enter] to continue")
  }
  cat("\n#################################", fill = TRUE)
  
  # Test fixneg argument
  cat("\nThe post-processing of negative predictions (should not give errors):", 
      fill = TRUE)
  dy_neg <- dy
  dy_neg[, 1:2] <- dy_neg[, 1:2] - 10
  m_t1_1_neg <- mgpr(datay = dy_neg, datax = dx, kernel = "matern32",
                 kernpar = list(sigma = 1, corlen = 10, errorvar = .1))
  cat("\nFixneg == FALSE: ", fill = TRUE)
  tryCatch({print(min(predict(m_t1_1_neg, kfold = 0, verbose = FALSE, 
                              fixneg = FALSE)))},
           error = function(e) print(e)) # NO CV
  cat("\nFixneg == TRUE: ", fill = TRUE)
  tryCatch({print(min(predict(m_t1_1_neg, kfold = 0, verbose = FALSE, 
                              fixneg = TRUE)))},
           error = function(e) print(e))  # NO CV
  cat("\nFixneg == FALSE (kfold): ", fill = TRUE)
  tryCatch({print(min(predict(m_t1_1_neg, kfold = 10, verbose = FALSE, 
                              fixneg = FALSE)))},
           error = function(e) print(e)) #KFOLD
  cat("\nFixneg == TRUE (kfold): ", fill = TRUE)
  tryCatch({print(min(predict(m_t1_1_neg, kfold = 10, verbose = FALSE, 
                              fixneg = TRUE)))},
           error = function(e) print(e)) # KFOLD
  cat("\nFixneg == FALSE (newdata): ", fill = TRUE)
  tryCatch({print(min(predict(m_t1_1_neg, newdatax = dx, verbose = FALSE, 
                              fixneg = FALSE)))},
           error = function(e) print(e)) #KFOLD
  cat("\nFixneg == TRUE (newdata): ", fill = TRUE)
  tryCatch({print(min(predict(m_t1_1_neg, newdatax = dx, verbose = FALSE, 
                              fixneg = TRUE)))},
           error = function(e) print(e)) # KFOLD
  
  if (manustep) {
    readline(prompt = "\nDoes it look good? Press [enter] to continue")
  }
  
  cat("\n#################################", fill = TRUE)
  
  # include also credible intervals
  cat("\nCheck manualy the credible intervals (should not give errors): ", 
      fill = TRUE)
  cat("\n", fill = TRUE)
  tryCatch({pred_lcv_ci <- predict(m_t1_1, verbose = F, kfold = 10, 
                                   credinter = 0.95)},
           error = function(e) print(e))
  tryCatch({pred_loocv_ci <- predict(m_t1_1, verbose = F, kfold = nrow(dx), 
                                     credinter = 0.95)},
           error = function(e) print(e))
  cat("Compare performance assessment: Kfold vs. LOOCV: ", fill = TRUE)
  print(validate(pred_lcv_ci$pred, dy))
  print(validate(pred_loocv_ci$pred, dy))
  
  cat("\nPredictions:", fill = TRUE)
  print(head(pred_lcv_ci$pred))
  print(head(pred_loocv_ci$pred))
  cat("\nCredible intervals:", fill = TRUE)
  print(head(pred_lcv_ci$credinter)) # lets see credible intervals
  print(head(pred_loocv_ci$credinter)) # lets see credible intervals
  
  cat("\nCheck negative credible intervals (fixneg = FALSE): ", fill = TRUE)
  cat("\nKFOLD mode: ", fill = TRUE)
  tryCatch(pred_lcv_ci_neg <- predict(m_t1_1_neg, verbose = F, kfold = 10, 
                                      credinter = 0.95))
  print(head(pred_lcv_ci_neg$pred))
  print(head(pred_lcv_ci_neg$credinter)) # lets see credible intervals
  
  cat("\nNew data mode: ", fill = TRUE)
  tryCatch(pred_lcv_ci_neg <- predict(m_t1_1_neg, verbose = F, newdatax = dx, 
                                      credinter = 0.95))
  print(head(pred_lcv_ci_neg$pred))
  cat("\n", fill = TRUE)
  print(head(pred_lcv_ci_neg$credinter)) # lets see credible intervals
  
  cat("\nCheck negative credible intervals (fixneg = TRUE): ", fill = TRUE)
  cat("\nKFOLD mode", fill = TRUE)
  tryCatch(pred_lcv_ci_negfix <- predict(m_t1_1_neg, verbose = F, kfold = 10, 
                                         credinter = 0.95, fixneg = TRUE))
  print(head(pred_lcv_ci_negfix$pred))
  cat("\n", fill = TRUE)
  print(head(pred_lcv_ci_negfix$credinter)) # lets see credible intervals
  
  cat("\nNew data mode", fill = TRUE)
  tryCatch(pred_lcv_ci_negfix <- predict(m_t1_1_neg, verbose = F, newdatax = dx, 
                                         credinter = 0.95, fixneg = TRUE))
  print(head(pred_lcv_ci_negfix$pred))
  cat("\n", fill = TRUE)
  print(head(pred_lcv_ci_negfix$credinter)) # lets see credible intervals
  
  if (manustep) {
    readline(prompt = "\nDoes it look good? Press [enter] to continue")
  }
  cat("\n#################################", fill = TRUE)
  
  # predictions with and without credible intervals should be equal
  cat(paste0("\nChecking that the LOOCV predictions are",  
             "similar with/without the credinter mode:"))
  if (mean(pred_m_t1_1_loocv == pred_loocv_ci$pred) == 1) {
    cat("\n", fill = TRUE)
    cat("Running....predictions are similar. Let's continue.", fill = TRUE)
  } else {
    stop("Error. Predictions are dissimilar.")
  }
  
  cat("\n#################################", fill = TRUE)
  # inputs as vectors; unimode
  cat("\nTesting a vector input for response (univariate):", fill = TRUE)
  cat("\n", fill = TRUE)
  tryCatch({m_t1_s <- mgpr( as.numeric(dy[, 1]), dx )},
           error = function(e) print(e))
  tryCatch({pred_s1 <- predict( m_t1_s, verbose = F)},
           error = function(e) print(e))
  tryCatch({summary(m_t1_s)},
           error = function(e) print(e))
  cat("\n", fill = TRUE)
  tryCatch({print(validate(pred_s1[,1], dy[,1]))},
    error = function(e) print(e))
  
  if (manustep) {
    readline(prompt = "Does it look good? Press [enter] to continue")
  }
  
  cat("\n#################################", fill = TRUE)
  # one response and one predictor: train & predict
  cat("\nTesting vector input for both response and predictor variable" ,
      fill = TRUE)
  tryCatch( {m_t1_s_s <- mgpr(as.numeric(dy[, 1]), dx[, 1])},
            error = function(e) print(e))
  tryCatch( {pred_s_s1 <- predict( m_t1_s_s, verbose=F)},
            error = function(e) print(e))
  tryCatch( {pred_s_s2 <- predict( m_t1_s_s, verbose=F, kfold = 10)},
            error = function(e) print(e))
  cat("\n", fill = TRUE)
  tryCatch({summary(m_t1_s_s)}, error = function(e) print(e))
  cat("\nPerformance assessment, train mode and kfold:", fill = TRUE)
  tryCatch({print(validate(pred_s_s1[, 1], dy[, 1]))},
           error = function(e) print(e))
    tryCatch({print(validate(pred_s_s2[, 1], dy[, 1]))},
             error = function(e) print(e))
  
  cat("\n#################################", fill = TRUE)
  cat("Testing with inputs that are not acceptable/reasonable:" ,
      fill = TRUE)
  cat("One training observation (should return an error notification):" ,fill = TRUE)
  tryCatch({m_dxna_oneobs <- mgpr(datay = dy[1, ], datax = dx[1, ], 
                                  kernel = "matern32",
                                  kernpar = list())},
                                  error = function(e) print(e))
  cat("\n", fill = TRUE)
  cat(paste0("Two training observations", 
          "(Do we wanna set a minimum here? How many? Based on the number of Xs?):"), 
      fill = TRUE)
  tryCatch({m_dxna_twoobs <- mgpr(datay = dy[1:2, ], datax = dx[1:2, ], 
                                  kernel = "matern32",
                                  kernpar = list())},
           error = function(e) print(e))
  tryCatch({summary(m_dxna_twoobs)},
           error = function(e) print(e))
  cat("\n", fill = TRUE)
  cat("\nTesting if NAs in the response data: " , fill = TRUE)
  dy_na <- dy
  dy_na[sample(1, nrow(dy), 3), 2] <- NA
  tryCatch({m_na_t <- mgpr(datay = dy_na, datax = dx, kernel = "matern32",
                           kernpar = list())},
           error = function(e) print(e)) 
  cat("\n", fill = TRUE)
  cat("Testing if NAs in the predictor variable data (should give errors) " ,
      fill = TRUE)
  dx_na <- dx
  dx_na[sample(1, nrow(dy), 3), 2] <- NA
  tryCatch({m_dxna_t <- mgpr(datay = dy, datax = dx_na, kernel = "matern32",
                             kernpar = list())},
           error = function(e) print(e))
  cat("\n", fill = TRUE)
  cat(paste0("Testing if var(x) == 0 in the predictor variable data", 
             " when using the newdata mode (should not return errors)"),
      fill = TRUE)
  dx_var0 <- dx
  dx_var0[, sample(1:ncol(dx_var0), 5)] <- 1
  cat("\nhead(newdata):", fill = TRUE)
  tryCatch({print(head(dx_var0))}, error = function(e) print(e))
  
  tryCatch({m_dxv0_t <- mgpr(datay = dy, datax = dx, kernel = "matern32",
                             kernpar = list())},
           error = function(e) print(e))
  tryCatch({pred_dx_var0 <- predict(m_dxv0_t, newdatax = dx_var0)},
           error = function(e) print(e))
  
  cat("\nhead(predictions):", fill = TRUE)
  tryCatch({head(pred_dx_var0)},
           error = function(e) print(e))
}


