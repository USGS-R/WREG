#' Weighing Matrix for ROI-WLS (WREG)
#'
#'@description
#'\code{Omega.WLS.ROImatchMatLab} calculates the weighting function for a 
#'weighted least-squares regression using regions-of-influence.  This is largely
#'legacy code to match WREG v. 1.05 idiosyncrasies.
#'
#' @param Y.all The dependent variable of interest at all sites in the network, with any
#' transformations already applied.
#' @param X.all The independent variables  at all sites in the network, with any 
#' transformations already applied.  Each row represents a site and each column 
#' represents a particular independent variable.  The rows must be in the same order as
#' the dependent variables in \code{Y}.
#' @param LP3 A dataframe containing the fitted Log-Pearson Type III standard
#' deviate, standard deviation and skew for all sites in the network.  The names of 
#' this data frame are \code{S}, \code{K} and \code{G}.  The order of the rows must 
#' be the same as \code{Y}.
#' @param RecordLengths.all \code{RecordLengths.all} may be a matrix whose rows and columns 
#' are in the same order as \code{Y}.  Each \code{(r,c)} element represents the 
#' length of concurrent record between sites \code{r} and \code{c}.  The diagonal 
#' elements therefore represent each site's full record length.  For \dQuote{WLS}, 
#' only the at-site record lengths are needed.  In this case, \code{RecordLengths} 
#' can be a vector or a matrix.  Here, this should only include sites in the current
#' region of influence.
#' @param NDX A vector listing the indices of the sites that comprise the region of 
#' influence.
#'
#'@details This is a legacy function that matches the idiosyncrasies of WREG v. 1.05.
#'This includes using all sites to implement the sigma regression, averaging across
#'all record lengths, using arbitrary record lengths to estimate weights and using a
#'new function for estimating the MSE of the basic OLS model.
#'
#'This function will become obsolete once all idiosyncrasies are assessed.
#'
#'@return \code{Omega.WLS.ROImatchMatLab} returns a list with three elements:
#'\item{Omega}{The estimated weighting matrix.}
#'\item{var.modelerror.0}{The estimated model error variance for a constant-value 
#'model.}
#'\item{var.modelerror.k}{The estimated model error variance for a k-variable model.}
#'
#'@export
Omega.WLS.ROImatchMatLab <- function(Y.all,X.all,LP3.all,RecordLengths.all,NDX) {
  # William Farmer, USGS, January 26, 2015
  #
  # WREG v 1.05 (MatLab) is inconsistent in which components are used to develop equations in RoI-WLS
  # Eq. 15 is implemented with all sites
  # The mean inverse record length is implemented with all sites
  
  ## Get the subset parameters
  Y<-Y.all[NDX] # Subset of dependent variables in region of influence.
  X<-X.all[NDX,] # Subset of independent variables in region of influence.
  LP3<-data.frame(LP3.all)[NDX,] # Subset of LP3 parameters for region of influence.
  if (is.matrix(RecordLengths.all)) {
    RecordLengths.all <- diag(RecordLengths.all)
  }
  RecordLengths <- RecordLengths.all[NDX]
  
  ## Initial OLS (basis for others)
  Omega <- diag(nrow(X)) # Temporary weighting matrix (identity) for initial OLS
  B_hat <- solve(t(X)%*%solve(Omega)%*%X)%*%t(X)%*%solve(Omega)%*%Y # OLS estimated coefficients
  Y_hat <- X%*%B_hat # OLS model estimates
  e <- Y-Y_hat # OLS model residuals
  MSE.OLS <- sum(e^2)/(nrow(X)-ncol(X)) # OLS mean square-error (k-variable)
  MSE.OLS <- sd(e)^2 # BUG: MatLab code uses this expression in RoI WLS.  See line 1450 of WREG v 1.05
  MSE.OLS.0 <- sum((Y-matrix(1,ncol=1,nrow=length(Y))%*%solve(t(matrix(1,ncol=1,nrow=length(Y)))%*%solve(Omega)%*%matrix(1,ncol=1,nrow=length(Y)))%*%t(matrix(1,ncol=1,nrow=length(Y)))%*%solve(Omega)%*%Y)^2)/(nrow(X)-1) # OLS mean square-error (constant model)
  
  ## Estimate k-variable model-error variance
  Omega <- diag(nrow(X.all)) # Temporary weighting matrix (identity) for k-variable sigma regression. Eq 15.
  B.SigReg <- solve(t(X.all)%*%solve(Omega)%*%X.all)%*%t(X.all)%*%solve(Omega)%*%LP3.all$S # OLS estimated coefficients for model of LP3 standard deviations
  # BUG: MatLab fits beta regression (Eq 15) using all sites.  See lines 1305-1316 and 1413-1419 of WREG v 1.05
  Yhat.SigReg <- X.all%*%B.SigReg # Estimates of sigma regression
  S_bar <- mean(Yhat.SigReg[NDX]) # Average sigma of limited set
  G_bar <- mean(LP3$G) # Average skew of limited-set LP3
  K_bar <- mean(LP3$K) # Average standard deviate of limited-set LP3
  c1 <- max(0,(1+K_bar^2*(1+0.75*G_bar^2)/2+K_bar*G_bar)*S_bar^2) # Coefficeint for calculation of the variance of at-site skew.  Eq 29.
  var.modelerror.k <- max(0,MSE.OLS-c1*mean(1/RecordLengths.all)) # k-variable model-error variance.  Eq 14.
  # BUG: MatLab code uses average of all inverse record lengths.  See lines 1453-1461.  
  
  ## Estimate 0-order model-error variance
  B.SigReg <- solve(t(matrix(1,ncol=1,nrow=nrow(X.all)))%*%solve(Omega)%*%matrix(1,ncol=1,nrow=nrow(X.all)))%*%t(matrix(1,ncol=1,nrow=nrow(X.all)))%*%solve(Omega)%*%LP3.all$S  # Temporary weighting matrix (identity) for constant-model sigma regression. Eq 15.
  Yhat.SigReg <- matrix(1,ncol=1,nrow=nrow(X.all))%*%B.SigReg # Estimates of sigma regression
  S_bar <- mean(Yhat.SigReg[NDX]) # Average sigma of limited set
  c1.0 <- max(0,(1+K_bar^2*(1+0.75*G_bar^2)/2+K_bar*G_bar)*S_bar^2) # Coefficeint for calculation of the variance of at-site skew.  Eq 29.
  var.modelerror.0 <- max(0,MSE.OLS.0-c1.0*mean(1/RecordLengths.all)) # Constant-model model-error variance.  Eq 14.
  
  ## Final weighting matrix
  Omega <- diag((var.modelerror.k+c1/RecordLengths.all[1:n])) # WLS weighting matrix.  Eq 12
  # BUG: MatLab code mis-shuffles the record lengths.  See lines 1457-1459 of WREG v. 1.05
  
  ## Output
  Output <- list(Omega=Omega,var.modelerror.0=var.modelerror.0,var.modelerror.k=var.modelerror.k)
  return(Output)
}
