#' Weighted-Multiple-Linear Regression Program (WREG)
#' 
#' @description
#' The \code{WREG.MLR} function executes the basic regression analysis that forms the 
#' foundation of the WREG program.
#' 
#' @param Y A numeric vector of the dependent variable of interest, with any transformations already 
#' applied.
#' @param X A numeric matrix of the independent variables in the regression, with any transformations 
#' already applied.  Each row represents a site and each column represents
#' a particular independe variable.  (If a leading constant is used, it should be 
#' included here as a leading column of ones.)  The rows must be in the same order as
#' the dependent variables in \code{Y}.
#' @param RecordLengths A numeric vector whose rows are in the same order as \code{Y} and represent the at-site record length.
#' @param LP3 A numeric matrix containing the fitted Log-Pearson Type III standard
#' deviate, standard deviation and skew for each site.  The columns of the matrix represent S, K, G, 
#' and an option regional skew value \code{GR} required by WREG.GLS with regSkew = TRUE. 
#' The order of the rows must be the same as \code{Y}.
#' @param x0 A vector containing the independent variables (as above) for a 
#' particular target site.  This variable is only used for ROI analysis.
#' @param Legacy A logical.  A value of \code{TRUE} forces the WREG program to behave 
#' identically to WREG v. 1.05, with BUGS and all.  It will force \code{TY=2} and
#' \code{DistMeth=1}.  For ROI regressions, it will also require a specific calculation 
#' for weighing matrices in \dQuote{WLS} (\code{\link{Omega.WLS.ROImatchMatLab}}), 
#' \dQuote{GLS}, and \dQuote{GLSskew} (see \code{\link{Omega.GLS.ROImatchMatLab}})
#' Further details are provided in \code{\link{WREG.RoI}}
#' 
#' @details
#' add details
#' 
#' @return All outputs are returned as part of a list.  The elements of the list depend 
#' on the type of regression performed.  The elements of the list may include:
#' \item{Coefs}{A data frame composed of four variables: (1) \code{Coefficient} contains 
#' the regression coefficeints estimated for the model, (2) \code{\sQuote{Standard Error}}
#' contains the standard errors of each regression coefficient, (3) \code{tStatistic} 
#' contains the Student's T-statistic of each regression coefficient and (4) 
#' \code{pValue} contains the significance probability of each regression coefficient.}
#' \item{ResLevInf}{A data frame composed of three variables for each site in the 
#' regression.  \code{Residual} contains the model residuals.  \code{Leverage} contains
#' the leverage of each site.  \code{Influence} contains the influence of each site.}
#' \item{LevLim}{The critical value of leverage.  See \code{\link{Leverage}}}
#' \item{InflLim}{The critical value of influence.  See \code{\link{Influence}}}
#' \item{LevInf.Sig}{A logical matrix indicating if the leverage (column 1) is significant 
#' and the influence (column 2) is significant for each site in the regression.}
#' \item{PerformanceMetrics}{A list of not more than ten elements.  All regression
#' types return the mean squared error of residuals (\code{MSE}), the coefficient of
#' determination (\code{R2}), the adjusted coefficient of determination (\code{R2_adj})
#' and the root mean squared error (\code{RMSE}, in percent).  \code{WLS}, \code{GLS} 
#' and \code{GLSskew} also return the pseudo coefficient of regression (\code{R2_pseudo}), 
#' the average variance of prediction (\code{AVP}), the standard error of prediction 
#' (\code{Sp}, in percent), a vector of the individual variances of prediction for 
#' each site (\code{VP.PredVar}), the model-error variance (\code{ModErrVar}) and the 
#' standardized model error variance (\code{StanModErr}, in percent).  Details on the 
#' appropriateness and applicability of performance metrics can be found in the WREG manual.}
#' \item{X}{The input predictors.}
#' \item{Y}{The input observations.}
#' \item{fitted.values}{A vector of model estimates from the regression model.}
#' \item{residuals}{A vector of model residuals.}
#' \item{Weighting}{The weighting matrix used to develop regression estimates.}
#' \item{Input}{A list of input parameters for error searching.  Right now it only
#' includes \code{Legacy} to document if a Legacy application of WREG v. 1.05 was
#' implemented.}
#' @import stats
#' 
#' @examples 
#' \dontrun Add example
#'@export

WREG.WLS <- function(Y,X,recordLengths,LP3,x0=NA,legacy=FALSE) {
  # William Farmer, USGS, January 05, 2015
  
  ## Add controls to meet legacy demands
  ##    NOTE: legacy forces program to return the same results as WREG v 1.05.
  if (legacy) { # If legacy is indicated, override custom inputs.
    TY <- 2 # WREG v1.05 does not read this input correctly.
    DistMeth <- 1 # WREG v1.05 uses "Nautical Mile" approximation
  }
  ## Determine if ROI is being applied
  if (is.na(sum(x0))) { # ROI regression is not used.
    ROI <- F
  } else { # ROI regression is used.
    ROI <- T
  }
  ## Just initial values for control.
  var.modelerror.k <- NA
  
  ### Correct recordLengths to be only at-site record lengths.
  if(is.matrix(recordLengths)) {
    recordLengths<-diag(recordLengths)
  }
  
  ###Rename LP3 columns
  ####Check if GR was supplied so that same LP3 file can be used for all functions
  if(ncol(LP3) == 3)
  {
  colnames(LP3) <- c("S","K","G")
  } else if(ncol(LP3) == 4)
  {
    colnames(LP3) <- c("S","K","G","GR")
  }
    
  #Convert X and Y from dataframes to matrices to work with matrix operations below
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  ### Initial OLS (basis for others)
  Omega <- diag(nrow(X)) # temporary weighting matrix (identity)
  B_hat <- solve(t(X)%*%solve(Omega)%*%X)%*%t(X)%*%solve(Omega)%*%Y # OLS estimated coefficients
  Y_hat <- X%*%B_hat # OLS model estimates
  e <- Y-Y_hat # OLS residuals
  MSE.OLS <- sum(e^2)/(nrow(X)-ncol(X)) # Mean square-error (k variable OLS)
  MSE.OLS.0 <- sum((Y-matrix(1,ncol=1,nrow=length(Y))%*%solve(t(matrix(1,ncol=1,nrow=length(Y)))%*%solve(Omega)%*%matrix(1,ncol=1,nrow=length(Y)))%*%t(matrix(1,ncol=1,nrow=length(Y)))%*%solve(Omega)%*%Y)^2)/(nrow(X)-1) # Mean square-error (0 variable, constant OLS)
  ### Estimate model-error variance
  B.SigReg <- solve(t(X)%*%solve(Omega)%*%X)%*%t(X)%*%solve(Omega)%*%LP3$S # OLS estimated coefficients for k-variable model of LP3 standard deviation. See Eq 15.
  Yhat.SigReg <- X%*%B.SigReg # Estimates of sigma regression
  S_bar <- mean(Yhat.SigReg) # average sigma of LP3
  G_bar <- mean(LP3$G) # average skew of LP3
  K_bar <- mean(LP3$K) # average standard deviate of LP3
  c1 <- max(0,(1+K_bar^2*(1+0.75*G_bar^2)/2+K_bar*G_bar)*S_bar^2) # Leading coefficient. Eq 13
  var.modelerror.k <- max(0,MSE.OLS-c1*mean(1/recordLengths)) # k-variable model-error variance. Eq 14.
  ### Estimate 0-order model-error variance
  B.SigReg <- solve(t(matrix(1,ncol=1,nrow=length(Y)))%*%solve(Omega)%*%matrix(1,ncol=1,nrow=length(Y)))%*%t(matrix(1,ncol=1,nrow=length(Y)))%*%solve(Omega)%*%LP3$S  # OLS estimated coefficients for constant model of LP3 standard deviation. See Eq 15.
  Yhat.SigReg <- matrix(1,ncol=1,nrow=length(Y))%*%B.SigReg # Estimates of sigma regression
  S_bar <- mean(Yhat.SigReg) # average sigma of LP3
  c1.0 <- max(0,(1+K_bar^2*(1+0.75*G_bar^2)/2+K_bar*G_bar)*S_bar^2) # Leading coefficient for constant-model model-error variance. See Eq 13
  var.modelerror.0 <- max(0,MSE.OLS.0-c1.0*mean(1/recordLengths)) # Constant-model model-error variance.  See Eq. 14.
  ### Final weighting matrix
  Omega <- diag((var.modelerror.k+c1/recordLengths)) # WLS weighting matrix.  Eq 12
  
  ## Basic regression calculations
  B_hat <- solve(t(X)%*%solve(Omega)%*%X)%*%t(X)%*%solve(Omega)%*%Y # Estimated regression coefficients. Eq 7 (9, 11, and 18)
  Y_hat <- X%*%B_hat # Model estimates. Eq 8
  e <- Y-Y_hat # Model residuals. Eq 30
  
  ## Performance metrics
  MSE <- sum(e^2)/(nrow(X)-ncol(X)) # Mean square-error. Eq 31
  SSR <- sum(e^2) # Residual sum of squares. Eq 36
  SST <- sum((Y-mean(Y))^2) # Total sum of squares. Eq 37
  R2 <- 1 - SSR/SST # Coefficient of determination. Eq 35
  R2_adj <- 1 - SSR*(nrow(X)-1)/SST/(nrow(X)-ncol(X)) # Adjusted coefficient of determination. Eq 38
  RMSE <-100*sqrt(exp(log(10)*log(10)*MSE)-1) # Root-mean-squared error, in percent. Eq 34
  PerfMet <- list(MSE=MSE,R2=R2,R2_adj=R2_adj,RMSE=RMSE) # Performance metrics for output (basic, for OLS)
  R2_pseudo <- 1 - var.modelerror.k/var.modelerror.0 # Pseudo coefficient of determination. Eq 39
  AVP <- var.modelerror.k + mean(diag(X%*%solve(t(X)%*%solve(Omega)%*%X)%*%t(X))) # Average varaince of prediction. Eq 32
  VP <- vector(length=length(Y)) # Empty vector for individual variances of prediction
  for (i in 1:length(VP)) { # Individual variances of prediction
    VP[i] <- var.modelerror.k + X[i,]%*%solve(t(X)%*%solve(Omega)%*%X)%*%X[i,] # Individual variance of prediction.  Based on Eq 32.
  }
  VP <- data.frame(VP); names(VP) <- 'PredVar' # Formating the VP vector for output
  Sp <- 100*sqrt(exp(log(10)*log(10)*AVP)-1) # Standard error of predictions. Eq 33.
  Se <-100*sqrt(exp(log(10)*log(10)*var.modelerror.k)-1) # Standard model error. Not noted in the manual, but included as output in WREG 1.05.  Based on Eq 33.
  PerfMet <- c(PerfMet,R2_pseudo=R2_pseudo,AVP=AVP,Sp=Sp,VP=VP,ModErrVar=var.modelerror.k,StanModErr=Se) # Performance metrics for output
  
  
  ## Leverage and influence statistics
  Lev <- Leverage(X=X,Omega=Omega,x0=x0,ROI=ROI) # Leverage subroutine
  Infl <- Influence(e=e,X=X,Omega=Omega,Beta=B_hat,ROI=ROI,Lev=Lev$Leverage) # Influence subroutine
  
  ## Significance of regression parameters
  B_var <- diag(solve(t(X)%*%solve(Omega)%*%X)) # Covariances of regression coefficients. Eq 46
  
  B_tval <- B_hat/sqrt(B_var) # T-value statistics of regression coefficients. Eq 45.
  B_pval <- 2*stats::pt(-abs(B_tval),df=(nrow(X)-ncol(X))) # Significnace of regression coefficients
  
  ## Create summary tables
  Coefs <- data.frame(cbind(B_hat,sqrt(B_var),B_tval,B_pval)) # Regression coefficient table for output
  names(Coefs) <- c('Coefficient','Standard Error','tStatistic','pValue')
  ResLevInf <- data.frame(cbind(e,Lev$Leverage,Infl$Influence)) # Residuals, leverage and influence of each varaible for output
  names(ResLevInf) <- c('Residual','Leverage','Influence')
  LevInf.Sig<-data.frame(cbind(Lev$Significant,Infl$Significant)) # Indication of significance for leverage and Influence for output
  names(LevInf.Sig) <- c('SignificantLeverage','SignificantInfluence')
  
  ## Handling output
  Output <- list(Coefs=Coefs,ResLevInf=ResLevInf,LevLim=Lev$Limit,InflLim=Infl$Limit,LevInf.Sig=cbind(Lev$Significant,Infl$Significant),PerformanceMetrics=PerfMet,X=X,Y=Y,fitted.values=Y_hat,residuals=e,Weighting=Omega,Inputs=list(legacy=legacy))
  if (ROI) { # Appended at-site estimates for ROI calculations
    Y_est <- x0%*%B_hat # ROI site estimate
    Output <- c(Output,Y.ROI=Y_est,x0.ROI=x0)
  }
  
  class(Output) <- 'WREG.WLS'
  return(Output)
}