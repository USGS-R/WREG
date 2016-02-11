#' Weighted-Multiple-Linear Regression Program (WREG)
#' 
#' @description
#' The \code{WREG.MLR} function executes the basic regression analysis that forms the 
#' foundation of the WREG program.
#' 
#' @param Y The dependent variable of interest, with any transformations already 
#' applied.
#' @param X The independent variables in the regression, with any transformations 
#' already applied.  Each row represents a site and each column represents
#' a particular independe variable.  (If a leading constant is used, it should be 
#' included here as a leading column of ones.)  The rows must be in the same order as
#' the dependent variables in \code{Y}.
#' @param x0 A vector containing the independent variables (as above) for a 
#' particular target site.  This variable is only used for ROI analysis.
#' @param Reg A string indicating which type of regression should be applied.
#' The options include: \dQuote{OLS} for ordinary least-squares regression,
#' \dQuote{WLS} for weighted least-squares regression, \dQuote{GLS} for generalized
#' least-squares regression, with no uncertainty in regional skew, \dQuote{GLSskew}
#' for generalized least-squares regression with uncertainty in regional skew, and 
#' \dQuote{CustomWeight} allowing the user to provide a custom weighting matrix.  
#' (In the case of \dQuote{GLSskew}, the uncertainty in regional skew must be provided
#'  as the mean squared error in regional skew.)
#' @param RecordLengths This input is required for \dQuote{WLS}, \dQuote{GLS} and 
#' \dQuote{GLSskew}.  For \dQuote{GLS} and \dQuote{GLSskew}, \code{RecordLengths} 
#' should be a matrix whose rows and columns are in the same order as \code{Y}.  Each 
#' \code{(r,c)} element represents the length of concurrent record between sites 
#' \code{r} and \code{c}.  The diagonal elements therefore represent each site's full
#' record length.  For \dQuote{WLS}, the only the at-site record lengths are needed.
#' In the case of \dQuote{WLS}, \code{RecordLengths} can be a vector or the matrix 
#' described for \dQuote{GLS} and \dQuote{GLSskew}.
#' @param LP3 A dataframe containing the fitted Log-Pearson Type III standard
#' deviate, standard deviation and skew for each site.  The names of this data frame are
#' \code{S}, \code{K} and \code{G}.  For \dQuote{GLSskew}, the regional skew value must 
#' also be provided in a variable called \code{GR}.  The order of the rows must be the same
#' as \code{Y}.
#' @param alpha A number, required only for \dQuote{GLS} and \dQuote{GLSskew}.  
#' \code{alpha} is a parameter used in the estimated cross-correlation between site
#' records.  See equation 20 in the WREG v. 1.05 manual.  The arbitrary, default value 
#' is 0.01.  The user should fit a different value as needed.
#' @param theta A number, required only for \dQuote{GLS} and \dQuote{GLSskew}.  
#' \code{theta} is a parameter used in the estimated cross-correlation between site
#' records.  See equation 20 in the WREG v. 1.05 manual.  The arbitrary, default value 
#' is 0.98.  The user should fit a different value as needed.
#' @param BasinChars A dataframe containing three variables: \code{StationID} is the 
#' numerical identifier (without a leading zero) of each site, \code{Lat} is the latitude
#' of the site, in decimal degrees, and \code{Long} is the longitude of the site, in decimal
#' degrees.  The sites must be presented in the same order as \code{Y}.  Required only for
#' \dQuote{GLS} and \dQuote{GLSskew}.
#' @param MSEGR A number. The mean squared error of the regional skew.  Required only for
#' \dQuote{GLSskew}.
#' @param TY A number.  The return period of the event being modeled.  Required only for 
#' \dQuote{GLSskew}.  The default value is \code{2}.  (See the \code{Legacy} details below.)
#' @param Peak A logical.  Indicates if the event being modeled is a peak flow event
#' or a low-flow event.  \code{TRUE} indicates a peak flow, while \code{FALSE} indicates
#' a low-flow event.
#' @param CustomWeight This allows the user to enter a custom weighting matrix.  It 
#' is included also to provide legacy code for WREG v. 1.05.  \code{CustomWeight} 
#' can either be a square matrix of weights with a length equal to \code{length(Y)}
#' or a list containing three elements.  The elements of the list include (1) 
#' \code{Omega} as the square weighting matrix, (2) \code{var.modelerror.k} as the 
#' estimated variance of the model errors from the k-variable model (\code{k=ncol(X)}),
#' and (3) \code{var.modelerror.0} as the variance of the model errors from a consant-only 
#' regression.  Required for \code{Reg=} \dQuote{CustomWeight}.
#' @param DistMeth Required for \dQuote{GLS} and \dQuote{GLSskew}.  A value of \code{1} 
#' indicates that the "Nautical Mile" approximation should be used to calculate inter-site
#' distances.  A value of \code{2} designates the Haversine approximation.  See 
#' \code{\link{Dist.WREG}}.  The default value is \code{2}.  (See the \code{Legacy} 
#' details below.)
#' @param Legacy A logical.  A value of \code{TRUE} forces the WREG program to behave 
#' identically to WREG v. 1.05, with BUGS and all.  It will force \code{TY=2} and
#' \code{DistMeth=1}.  For ROI regressions, it will also require a specific calculation 
#' for weighing matrices in \dQuote{WLS} (\code{\link{Omega.WLS.ROImatchMatLab}}), 
#' \dQuote{GLS}, and \dQuote{GLSskew} (see \code{\link{Omega.GLS.ROImatchMatLab}})
#' Further details are provided in \code{\link{WREG.RoI}}
#' 
#' @details
#' This is the main function for WREG.  It builds weighing matrices, applies the specified
#' regression and computes the relevant performance metrics. Additional information on this 
#' program can be found in the manual of WREG v. 1.0.
#' 
#' The logical handle \code{Legacy} has been included to test that this program returns the 
#' same results as WREG v. 1.05.  In the development of this code, some idiosyncrasies of 
#' the MatLab code for WREG v. 1.05 became apparent.  Setting \code{Legacy} equal to 
#' \code{TRUE} forces the code to use the same idiosycrasies as WREG v. 1.05.  Some of 
#' these idosyncrasies may be bugs in the code.  Further analysis is needed.  For 
#' information on the specific idiosyncrasies, see the notes for the \code{Legacy} 
#' input and the links to other functions in this package.
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
#'data("baseData",package="WREG")
#'names(baseData$Dependent[2:4])
#'names(baseData$Independent[9:14])
#'exY <- log(baseData$Dependent$Q1.)
#'X1 <- log(baseData$Independent$DRNAREA)
#'X2 <- log(baseData$Independent$PRECIP)
#'X0 <- rep(1,length(X1))
#'exX <- cbind(X0,X1,X2)
#'Ex.OLS <- WREG.MLR(Y=exY,X=exX,Reg='OLS')
#'@export
WREG.MLR <- function(Y,X,x0=NA,Reg=c('OLS','WLS','GLS','GLSskew','CustomWeight'),
  RecordLengths=NA,LP3=NA,alpha=0.01,theta=0.98,BasinChars=NA,
  MSEGR=NA,TY=2,Peak=T,CustomWeight=NA,DistMeth=2,Legacy=FALSE) {
  # William Farmer, USGS, January 05, 2015
  
  ## Add controls to meet Legacy demands
  ##    NOTE: Legacy forces program to return the same results as WREG v 1.05.
  if (Legacy) { # If legacy is indicated, override custom inputs.
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
  CustomModelError <- FALSE
  
  ## Weighting matrix
  if (Reg=='OLS') { # Ordinary Least Squares
    Omega <- diag(nrow(X)) # weighting matrix
  } else if (Reg=='WLS') { # for WLS
    ### Correct RecordLengths to be only at-site record lengths.
    if(is.matrix(RecordLengths)) {
      RecordLengths<-diag(RecordLengths)
    }
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
    var.modelerror.k <- max(0,MSE.OLS-c1*mean(1/RecordLengths)) # k-variable model-error variance. Eq 14.
    ### Estimate 0-order model-error variance
    B.SigReg <- solve(t(matrix(1,ncol=1,nrow=length(Y)))%*%solve(Omega)%*%matrix(1,ncol=1,nrow=length(Y)))%*%t(matrix(1,ncol=1,nrow=length(Y)))%*%solve(Omega)%*%LP3$S  # OLS estimated coefficients for constant model of LP3 standard deviation. See Eq 15.
    Yhat.SigReg <- matrix(1,ncol=1,nrow=length(Y))%*%B.SigReg # Estimates of sigma regression
    S_bar <- mean(Yhat.SigReg) # average sigma of LP3
    c1.0 <- max(0,(1+K_bar^2*(1+0.75*G_bar^2)/2+K_bar*G_bar)*S_bar^2) # Leading coefficient for constant-model model-error variance. See Eq 13
    var.modelerror.0 <- max(0,MSE.OLS.0-c1.0*mean(1/RecordLengths)) # Constant-model model-error variance.  See Eq. 14.
    ### Final weighting matrix
    Omega <- diag((var.modelerror.k+c1/RecordLengths)) # WLS weighting matrix.  Eq 12
  } else if (is.element(Reg,c('GLS','GLSskew'))) { # for GLS and GLS skew
    if (Reg=='GLS') {MSEGR<-NA}
    ### Estimate k-variable model-error variance and final weighting matrix
    GLS.Weighting <- Omega.GLS(alpha=alpha,theta=theta,Independent=BasinChars,X=X,Y=Y,RecordLengths=RecordLengths,LP3=LP3,MSEGR=MSEGR,TY=TY,Peak=Peak,DistMeth=DistMeth) # Subroutine to calculate GLS weighting by iteration
    var.modelerror.k <- GLS.Weighting$GSQ # k-variable model-error variance
    Omega <- GLS.Weighting$Omega # GLS/GLSskew weighting matrix
    ### Estimate 0-order model-error variance
    X.0 <- matrix(1,ncol=1,nrow=length(Y)) # Constant predictor
    GLS.Weighting <- Omega.GLS(alpha=alpha,theta=theta,Independent=BasinChars,X=X.0,Y=Y,RecordLengths=RecordLengths,LP3=LP3,MSEGR=MSEGR,TY=TY,Peak=Peak,DistMeth=DistMeth) # Subroutine to calculate GLS weighting by iteration
    var.modelerror.0 <- GLS.Weighting$GSQ # constant-model model-error variance
  } else if (Reg=='CustomWeight') { # Allows user to specify particular weighting scheme.  Useful for legacy code.
    if (is.list(CustomWeight)) { # Omega from legacy code; Also contains information on model-error variances.
      Omega <- CustomWeight$Omega # Custom weighting
      var.modelerror.k <- CustomWeight$var.modelerror.k # Custom k-variable model-error variance
      var.modelerror.0 <- CustomWeight$var.modelerror.0 # Custom constant-model model-error variance
      CustomModelError=TRUE # Logical to note that the user has provided information on the model-error variance
    } else { # User-defined weighting, with no informaiton on model-error variance.
      Omega <- CustomWeight # Custom weighting
      var.modelerror.k <- NA # NULL custom k-variable model-error variance to control for errors.
      var.modelerror.0 <- NA # NULL custom constant-model model-error variance to control for errors.
    }
  }
  
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
  if (is.element(Reg,c('WLS','GLS','GLSskew'))||CustomModelError==T) { # non-OLS requires additional performance metrics
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
  }
  
  ## Leverage and influence statistics
  Lev <- Leverage(X=X,Omega=Omega,x0=x0,ROI=ROI) # Leverage subroutine
  Infl <- Influence(e=e,X=X,Omega=Omega,Beta=B_hat,ROI=ROI,Lev=Lev$Leverage) # Influence subroutine
  
  ## Significance of regression parameters
  B_var <- diag(solve(t(X)%*%solve(Omega)%*%X)) # Covariances of regression coefficients. Eq 46
  if (Reg=='OLS') { # Eq 46 in Manual for v1.05 is incoorect.  As reflected in code for v1.05 and independent verification, formula is altered for OLS.
    B_var <- MSE*B_var # Altered Eq 46 for OLS
  }
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
  Output <- list(Coefs=Coefs,ResLevInf=ResLevInf,LevLim=Lev$Limit,InflLim=Infl$Limit,LevInf.Sig=cbind(Lev$Significant,Infl$Significant),PerformanceMetrics=PerfMet,X=X,Y=Y,fitted.values=Y_hat,residuals=e,Weighting=Omega,Inputs=list(Legacy=Legacy))
  if (ROI) { # Appended at-site estimates for ROI calculations
    Y_est <- x0%*%B_hat # ROI site estimate
    Output <- c(Output,Y.ROI=Y_est,x0.ROI=x0)
  }
  
  if (Reg=='OLS') {
    class(Output) <- 'WREG.OLS'
  } else if (Reg=='WLS') {
    class(Output) <- 'WREG.WLS'
  } else if (Reg=='GLS') {
    class(Output) <- 'WREG.GLS'
  } else if (Reg=='GLSskew') {
    class(Output) <- 'WREG.GLSs'
  }
  return(Output)
}