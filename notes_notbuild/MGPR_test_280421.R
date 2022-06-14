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

#
# Set path
#
require(devtools) # Install if not alreary installed
polku <- "C:\\Users\\janr\\OneDrive - Norsk Institutt for Bioøkonomi\\MGPR\\MGPR_dev\\Rpaketti\\MGPR_rlib"
install_local(path = paste0(polku),
              force = TRUE)

require(MGPR)

#
# Benchmark helper funtion
#
benchmark <- function( obj )
{
  time1 = Sys.time()
  obj
  time2 = Sys.time()
  #out = round(as.numeric(paste(time2-time1)),3)
  out = round(difftime(time2,as.POSIXct(time1),units="secs"),3) #JR added: this returns always seconds

  cat( "Execution time", out, "seconds", "\n" )
}

#
# Performance assessment function. Inputs are data.frames
# having corresponding columns with matching column names.
# Returns absolute and relative rmse and bias inside of a
# data.frame.
#
validate <- function( pred, obs )
{
  if ( is.numeric(pred) & is.numeric(obs) )
  {
     pred = as.data.frame( pred )
     obs = as.data.frame( obs )
     colnames(pred) = colnames(obs) = "unknown"
  }

  if ( !is.data.frame(pred) | !is.data.frame(obs) )
  {
    print("ERROR: pred and obs must be of type 'data.frame'.")
    return( NULL )
  }
  if ( mean(colnames(pred)==colnames(obs)) != 1 )
  {
    print("ERROR: colnames in pred and obs in validate do not match.")
    return( NULL )
  }
  m = matrix( 0, nrow=4, ncol=ncol(pred) )

  for( i in 1:ncol(pred))
  {
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
# Test case 1: 100 obs, 3 ys and 25 xs.
#
data(mgprdata)

dx = mgprdata[, 6:29]
dy = mgprdata[, 1:5]


# train & predict
benchmark( m_t1 <<- mgpr( dy, dx) ) # Error
benchmark( m_t1 <<- mgpr(datay = dy, datax = dx, kernel = "matern3")) # error
benchmark( m_t1 <<- mgpr(datay = dy, datax = dx, kernel = "matern32")) # Works
summary( m_t1 ) # default kernel hpars

# Define kern hpars manually
m_t2 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
             kernpar = list()) # Throws an error
m_t2 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
             kernpar = list(corlen = 1)) # Throws an error since both hpars needed
m_t2 <- mgpr(datay = dy, datax = dx, kernel = "matern32",
             kernpar = list(corlen = 5, sigma = 1)) # works
summary( m_t2 ) # manually defined kernel hpars

benchmark( pred_t2 <<- predict( m_t1, verbose = T, kfold = 10,
                                covout = T, credinter = .95)) 

names(pred_t1)
pred_t1 <- predict( m_t1, verbose = F, covout = F)
benchmark( pred_t2 <<- predict( m_t2, verbose = F ) )
validate( pred_t1, dy )
validate( pred_t2, dy )

# Try Radial basis kernel
m_t3 <- mgpr(datay = dy, datax = dx, kernel = "rbf", kernpar = list()) # Throws an error
m_t3 <- mgpr(datay = dy, datax = dx, kernel = "rbf", kernpar = list(corlen = 1)) # Throws an error since both hpars needed
m_t3_1 <- mgpr(datay = dy, datax = dx, kernel = "rbf", kernpar = list(corlen = 3, sigma = 1)) # works; HOX! corlen 2; overfits
m_t3_2 <- mgpr(datay = dy, datax = dx, kernel = "rbf", kernpar = list(corlen = 10, sigma = 1)) # works; HOX! corlen 2; overfits

summary( m_t3_1)
summary( m_t3_2) # manually defined kernel hpars

benchmark( pred_t3_1 <<- predict( m_t3_1, verbose = F ) )
benchmark( pred_t3_2 <<- predict( m_t3_2, verbose = F ) )
validate( pred_t3_1, dy ) # HOX! corlen 2 with RBF; "overfits"
validate( pred_t3_2, dy ) # HOX! corlen 10 with RBF; much "realistic" in terms of overfitting


# kfold (kfold = 0 is default; kfold = 1 is not allowed; otherwise CVs)
pred_kfold <- predict( m_t1, kfold = 10, verbose = FALSE)
validate( pred_kfold, dy )
pred_kfold_2 <- predict( m_t2, kfold = 10, verbose = FALSE)
validate( pred_kfold_2, dy )
pred_kfold_3 <- predict( m_t3_2, kfold = 10, verbose = TRUE)
validate( pred_kfold_2, dy )

# Test fixneg argument
min(predict( m_t2, kfold = 0, verbose = FALSE, fixneg = FALSE)) # NO CV
min(predict( m_t2, kfold = 0, verbose = FALSE, fixneg = TRUE))  # NO CV
min(predict( m_t2, kfold = 10, verbose = FALSE, fixneg = TRUE)) # KFOLD
min(predict( m_t2, kfold = 10, verbose = FALSE, fixneg = FALSE)) #KFOLD

# loocv
pred_lcv <<- predict( m_t1, kfold = 100, verbose = FALSE)
validate( pred_lcv, dy )

# include also credible intervals
benchmark( pred_lcv_ci <<- predict( m_t1, verbose = T, kfold = 100, credinter=0.95 ) )
validate( pred_lcv_ci$pred, dy )
print(head(pred_lcv_ci$credinter)) # lets see credible intervals

# predictions with and without credible intervals should be equal
mean(pred_lcv == pred_lcv_ci$pred) == 1

# kfold with credible intervals
benchmark( pred_kfold_ci <<- mgpr_kfold( m_t1, 10, verbose=F, credinter=0.95 ) )
validate( pred_kfold_ci[,1:3], dy )

# loocv with credible intervals
benchmark( pred_lcv_ci <<- mgpr_kfold( m_t1, nrow(dy), verbose=F, credinter=0.95 ) )
validate( pred_lcv_ci[,1:3], dy )

# inputs as vectors
benchmark( m_t1_s <<- mgpr_train( as.numeric(dy[,1]), dx ) )
benchmark( pred_s1 <<- predict( m_t1_s, verbose=F ) )
summary( m_t1_s )
validate( pred_s1[,1], dy[,1] )

# one response and one predictor: train & predict
benchmark( m_t1_s_s <<- mgpr_train( as.numeric(dy[,1]), dx[,1], normPred=T ) )
benchmark( pred_s_s <<- predict( m_t1_s_s, verbose=F ) )
summary( m_t1_s_s )
validate( pred_s_s[,1], dy[,1] )

# one response and one predictor: kfold
benchmark( m_t1_s_s_kfold <<- mgpr_kfold( m_t1_s_s, 10, verbose=F ) )
validate( m_t1_s_s_kfold[,1], dy[,1] )

# one response and one predictor: loocv
benchmark( m_t1_s_s_loo <<- mgpr_kfold( m_t1_s_s, length(dx[,1]), verbose=F ) )
validate( m_t1_s_s_loo[,1], dy[,1] )

# predict to newdata (original t1 data is splitted)
# set.seed( 133 )
# samp = sample( seq(1,100,1), 30 )
# new_dx = dx[ samp, ]
# y_new = dy[ samp, ]
# old_dx = dx[ -samp, ]
# y_old = dy[ -samp, ]
# 
# benchmark( m_t1_2 <<- mgpr_train( y_old, old_dx, normPred=T ) )
# benchmark( pred_t1_2 <<- predict( m_t1_2, verbose=F ) )
# summary( m_t1_2 )
# 
# benchmark( pred_tonew <<- predict( m_t1_2, verbose=F, newData=new_dx) )
# validate( pred_tonew, y_new  )
# 
# benchmark( pred_tonew_ci <<- predict( m_t1_2, verbose=F, newData=new_dx, normPred=T, credinter=0.95) )
# validate( pred_tonew_ci[,1:3], y_new  )
# head( pred_tonew_ci[,4:6] )
# head( pred_tonew_ci[,7:9] )
# 
# # lets try if predictor var() == 0
# dx2 = dx; dx2[,c(10,12)] = rep( 1, nrow(dx) )
# benchmark( m_t_var <<- mgpr_train( dy, dx2, normPred=T ) ) # should fail
# 
# # what do the function return if we add total volume to the response set?
# dy$V_tot = rowSums( dy )
# benchmark( m_t_totV <<- mgpr_train( dy, dx, normPred=T ) )
# benchmark( pred_totV <<- predict( m_t_totV, verbose=F ) )
# summary( m_t_totV )
# validate( pred_totV, dy )
# pred_lcv_vtot <<- mgpr_kfold( m_t_totV, nrow(dy), verbose=F )
# validate( pred_lcv_vtot, dy )

#
# Test case 2: 493 obs, 15 ys and 77 xs. LOOCV run is exactly the
# same as in Varvia et al. (2017)
#
# 
# dx = read.table( paste( test_path, "testdata493_x77.txt", sep="" ), header=FALSE )
# dy = read.table( paste( test_path, "testdata493_y15.txt", sep="" ), header=FALSE )
# 
# benchmark( m_t2 <<- mgpr_train( dy, dx, normPred=T ) )
# benchmark( pred <<- predict( m_t2, verbose=F ) ) # takes app. 1 minute
# summary( m_t2 )
# 
# # loocv results should match with the Varvia et al. 2017 (takes app. 100 minutes)
# benchmark( pred_lcv <<- mgpr_kfold( m_t2, k=nrow(dy), verbose=F, credinter=0.95 ) )
# validate( pred_lcv[,1:15], dx )
