# Functions for WREG in R
# William Farmer

# Systemic changes
# 01/29/2015, WHF: Revised formulation of "RecordLengths" so that 
#                 "FileDir" is no longer calling outside files.
#                 Also revised Omega.WLS.ROImatchMatLab to have fewer inputs.

#' Leverage Statistics (WREG)
#'
#' @description
#' The \code{Leverage} function calculates the leverage statistics for each observation 
#' in the regression.
#' 
#' @param X contains independent variables.  A leading constant of one is included if 
#' a constant is included in the regression.  Each row represents a unique observation.
#' @param Omega is the weighting matrix used for regression fitting
#' @param Ch allows the user to specify a custom criteria for testing.  If not specified,
#'  the default is 4 for region-of-influence regression and 2 elsewise.
#' @param x0 is used only in the case of region-of-influence regression.  It contains the
#'  independent variables of a specific site against which to calculate the leverage.
#' @param ROI is a logical vector specifying if the regression is or is not 
#' region-of-influence regression
#' 
#' @details
#' Leverage is a measure of how far way the independent variables of an observation 
#' are from the other observations in the regression set.  A leverage is considered 
#' significant if the absolute value of the leverage is greater than the critical 
#' value.  These calculations are based on equations 40, 41 and 42 of the WREG v. 1.0 manual.
#' 
#' @return The function returns a list as output.  This list contains:
#' \item{Leverage}{A vector indicating the leverage of each observation.  
#' This is a vector whose length equals the number of rows in X.}
#' \item{Limit}{A number indicating the critical value of leverage for this data set.}
#' \item{Significant}{A logical vector the same size as \code{Leverage}. 
#'  It indicates if the leverage is significant for each observation.}
#'@export
Leverage <- function(X,Omega,Ch=NA,x0=NA,ROI=c(TRUE,FALSE)) {
  # William Farmer, USGS, January 02, 2015
  # 01/27/2015, WHF: Added ABS on limits
  
  if (ROI==F) { # Not region-of-influence regression
    h <- diag(X%*%solve(t(X)%*%solve(Omega)%*%X)%*%t(X)%*%solve(Omega)) # Leverage, Eq 40
    if (is.na(Ch)) {Ch=2} # Default for non-ROI 
  } else if (ROI==T) { # Region-of-influence regression
    h <- t(x0%*%solve(t(X)%*%solve(Omega)%*%X)%*%t(X)%*%solve(Omega)) # Leverage, Eq 41
    if (is.na(Ch)) {Ch=4} # Default for ROI 
  }
  h_limit <- Ch*mean(h) # Critical value of leverage, Eq 42
  sig <- abs(h)>abs(h_limit) # Significance test. Logical: Leverage is signifianct.
  Output <- list(Leverage=h,Limit=h_limit,Significant=sig)
  return(Output)
}

#' Influence Statistics (Cook's D) (WREG)
#' 
#' @description
#' The \code{Influence} function calculates the influence statistics (Cooks-D) for any 
#' regrerssion.
#'
#' @param e contains the model residuals.
#' @param X contains independent variables.  A leading constant of one is included if 
#' a constant is included in the regression.  Each row represents a unique obsevation.
#' @param Omega is the weighting matrix used for regression fitting.
#' @param Beta contains the fitted model coefficients.
#' @param ROI is a logical indicating if this is a region-of-influence regression.
#' @param Lev is a vector with the same length as \code{e} and includes the leverage of 
#' each observation.  
#' This input is required for any region-of-influence regression.
#' 
#' @details
#' Influence is a measure of the impact each observation has on the estimated 
#' regression coefficients.  The calculation is based on equation 43 of the WREG 
#' v. 1.0 manual.  The critical value of influence is calculated using equation 44.  
#' An influence is considered significant if the absolute value of the influence 
#' is greater than the critical value.
#' 
#' For region-of-influence regressions, the influence calculation is weighted 
#' by the leverage of that observation on the target site and the overall leverage 
#' of the observation.  This is a departure from the WREG v. 1.0, but reflects the 
#' WREG v. 1.05 code.
#' 
#' @return The function returns as list as output.  The list contains:
#'\item{Influence}{A vector containing the influence (Cook's D) of each observation 
#'on the estimated regression coefficients.}
#'\item{Limit}{The critical influence value for this dataset.}
#'\item{Significant}{A logical vector the same size as \code{Influence}. 
#'  It indicates if the influence is significant for each observation.}
#'@export
Influence <- function(e,X,Omega,Beta,ROI=FALSE,Lev=NA) {
  # William Farmer, USGS, January 02, 2015
  # 01/27/2015, WHF: Added abs on limit and ROI option to match MatLab WREG v 1.05, RoIMetrics, Lines 39-52
  
  L <- X%*%solve(t(X)%*%solve(Omega)%*%X)%*%t(X) # Basic leverage calculation
  if (ROI) {
    # Though not noted in the manual, GLS-ROI uses a different calculation and critical value for influence.
    H <- X%*%solve(t(X)%*%solve(Omega)%*%X)%*%t(X)%*%solve(Omega) # Leverage matrix; need for ROI version of leverage
    D <- (e^2)*diag(L)*abs(Lev/diag(H))/length(Beta)/(diag(Omega)-diag(L))^2 # ROI influence
    D_limit <- 4*2/nrow(X) # ROI critical influence
  } else {
    # Influence formula noted in manual.  Only used for non-ROI.
    D <- (e^2)*diag(L)/length(Beta)/(diag(Omega)-diag(L))^2 # Influence, Eq 43
    D_limit <- 4/nrow(X) # Critical influence, Eq 44
  }
  sig <- abs(D)>abs(D_limit) # Significance test. Logical: Influence is signifianct.
  Output <- list(Influence=D,Limit=D_limit,Significant=sig)
  return(Output)
}

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
  B_pval <- 2*pt(-abs(B_tval),df=(nrow(X)-ncol(X))) # Significnace of regression coefficients
  
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
  return(Output)
}

#' Calculate weighting matrix for GLS regression. (WREG)
#'
#'@description 
#'THe \code{Omega.GLS} function calculates the weighting matrix required for generalized
#'least-squares regression, without or without uncertainty in the regional skew.
#'
#' @param alpha A number, required only for \dQuote{GLS} and \dQuote{GLSskew}.  
#' \code{alpha} is a parameter used in the estimated cross-correlation between site
#' records.  See equation 20 in the WREG v. 1.05 manual.  The arbitrary, default value 
#' is 0.01.  The user should fit a different value as needed.
#' @param theta A number, required only for \dQuote{GLS} and \dQuote{GLSskew}.  
#' \code{theta} is a parameter used in the estimated cross-correlation between site
#' records.  See equation 20 in the WREG v. 1.05 manual.  The arbitrary, default value 
#' is 0.98.  The user should fit a different value as needed.
#' @param Independent A dataframe containing three variables: \code{StationID} is the 
#' numerical identifier (without a leading zero) of each site, \code{Lat} is the latitude
#' of the site, in decimal degrees, and \code{Long} is the longitude of the site, in decimal
#' degrees.  The sites must be presented in the same order as \code{Y}.  Required only for
#' \dQuote{GLS} and \dQuote{GLSskew}.
#' @param X The independent variables in the regression, with any transformations 
#' already applied.  Each row represents a site and each column represents
#' a particular independe variable.  (If a leading constant is used, it should be 
#' included here as a leading column of ones.)  The rows must be in the same order as
#' the dependent variables in \code{Y}.
#' @param Y The dependent variable of interest, with any transformations already 
#' applied.
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
#' @param MSEGR A number. The mean squared error of the regional skew.  Required only for
#' \dQuote{GLSskew}.
#' @param TY A number.  The return period of the event being modeled.  Required only for 
#' \dQuote{GLSskew}.  The default value is \code{2}.  (See the \code{Legacy} details below.)
#' @param Peak A logical.  Indicates if the event being modeled is a peak flow event
#' or a low-flow event.  \code{TRUE} indicates a peak flow, while \code{FALSE} indicates
#' a low-flow event.
#' @param DistMeth Required for \dQuote{GLS} and \dQuote{GLSskew}.  A value of \code{1} 
#' indicates that the "Nautical Mile" approximation should be used to calculate inter-site
#' distances.  A value of \code{2} designates the Haversine approximation.  See 
#' \code{\link{Dist.WREG}}.  The default value is \code{2}.  (See the \code{Legacy} 
#' details below.)
#'
#'@details This function is largely a subroutine for \code{\link{WREG.MLR}} when applying
#'generalized least-squares regression (\dQuote{GLS} and \dQuote{GLSskew}).
#'
#'The weighting matrix is calculated by iteration, as noted in the manual to 
#'WREG v. 1.0.  As currently implemented the initial estimate of model error variance,
#'\code{GSQ}, is taken to range from \code{0} and \code{2*var(Y)}.  This interval is 
#'broken into 30 equally spaced intervals.  The weighting matrix is calculated for each
#'interval endpoint and the deviation from equation 21 in the WREG v. 1.0 manual is
#'recorded.  The progam then search for the interval over which the deviatioon changes 
#'sign.  This interval is then split into 30 finer intervals and the process is repeated.
#'The ine interval with the smallest positive deviation is selected as the best estimate.
#'
#'@return This function returns a list with two elements:
#'\item{GSQ}{The estimated model error variance.}
#'\item{Omega}{The estimated weighting matrix.  A square matrix.}
#'@export
Omega.GLS <- function(alpha=0.01,theta=0.98,Independent,X,Y,RecordLengths,
  LP3,MSEGR=NA,TY=2,Peak=T,DistMeth=2) {
  # William Farmer, January 22, 2015
  
  ## Determining if skew adjustment is requested
  SkewAdj<-F # default: no skew adjustment
  if (!is.na(MSEGR)) {
    SkewAdj<-T # If user provides a mean squared-error of regional skew, 
    #               then use skew adjustment.
  }
  
  ## Create distance matrix and concurrent record lengths
  ##    (For skew-adjusted GLS: Also calculates the LP3 partial derivatives, mean squared-errors of at-site skew and the variance of at-site skew.)
  Dists <- matrix(NA,ncol=length(Y),nrow=length(Y)) # Empty matrix for intersite distances
  M <- RecordLengths # Just renaming input for ease.  (Should probably correct later...)
  if (SkewAdj) { # Make empty vectors for GLS-skew
    dKdG <- vector(length=length(Y)) # Empty vector for LP3 partial derivatives
    MSEg <- vector(length=length(Y)) # Empty vector for mean squared-error of at-site skew
    Varg <- vector(length=length(Y)) # Empty vector for variance of at-site skew
    ### Convert return period into probability
    if (Peak) { # if a peak flow is being estimated
      Zp <- -qnorm(1/TY) 
    } else { # if a low flow is being estimated
      Zp <- qnorm(1/TY)
    }
  }
  for (i in 1:length(Y)) {
    for (j in 1:length(Y)) {
      if (i!=j) {
        ### Calculate intersite distance via subroutine.
        ###     DistMeth==1 applies the 'Nautical Mile' approximation from WREG v 1.05
        ###     DistMeth==2 applies Haversine approximation
        Dists[i,j] <- Dist.WREG(Lat1=Independent$Lat[i],Long1=Independent$Long[i],Lat2=Independent$Lat[j],Long2=Independent$Long[j],method=DistMeth) # Intersite distance, miles
      }
    }
    if (SkewAdj) { # Additional calculations for skew-adjusted GLS
      ### LP3 Partial Derivative
      dKdG[i] <- (Zp^2-1)/6+LP3$G[i]*(Zp^3-6*Zp)/54-LP3$G[i]^2*(Zp^2-1)/72+Zp*LP3$G[i]^3/324+5*LP3$G[i]^4/23328 # LP3 partial derivative. Eq 23.
      ### Variance of at-site skew
      a <- -17.75/M[i,i]^2+50.06/M[i,i]^3 # Coefficeint for calculation of the variance of at-site skew.  Eq 27.
      b1 <- 3.92/M[i,i]^0.3-31.1/M[i,i]^0.6+34.86/M[i,i]^0.9  # Coefficeint for calculation of the variance of at-site skew.  Eq 28.
      c1 <- -7.31/M[i,i]^0.59+45.9/M[i,i]^1.18-86.5/M[i,i]^1.77  # Coefficeint for calculation of the variance of at-site skew.  Eq 29.
      Varg[i] <- (6/M[i,i]+a)*(1+LP3$GR[i]^2*(9/6+b1)+LP3$GR[i]^4*(15/48+c1)) # Variance of at-site skew. Eq 26.
      ### Mean squared-error of at-site skew
      ###   These equations are not included or described in v1.05 manual.  They are similar to Eq 28 and 29.
      ###   These equations, which appear in the v1.05 code, are documented in Griffis and Stedinger (2009); Eq 3, 4, 5 and 6.
      b <- 3.93/M[i,i]^0.3-30.97/M[i,i]^0.6+37.1/M[i,i]^0.9 # Coefficient for the calculation of mean squared-error of at-site skew.
      c <- -6.16/M[i,i]^0.56+36.83/M[i,i]^1.12-66.9/M[i,i]^1.68 # Coefficient for the calculation of mean squared-error of at-site skew.
      MSEg[i] <- (6/M[i,i]+a)*(1+LP3$G[i]^2*(9/6+b)+LP3$G[i]^4*(15/48+c)) # Mean squared-error of at-site skew.
    }
  }
  Rhos <- theta^(Dists/(alpha*Dists+1)) # Estimated intersite correlation. Eq 20.
  
  if (SkewAdj) { # if skew-adjusted GLS
    ### Calculate skew weights.
    Wg <- MSEGR/(MSEg+MSEGR) # skew weight. Eq 17.
    ### Calculate covariances between at-site skews
    Covgg <- matrix(NA,ncol=length(Y),nrow=length(Y)) # Empty matrix for covariances between at-site skews.
    for (i in 1:length(Y)) {
      for (j in 1:length(Y)) {
        if (i!=j) {
          Covgg[i,j] <- M[i,j]*sign(Rhos[i,j])*abs(Rhos[i,j])^3*sqrt(Varg[i]*Varg[j])/sqrt((M[i,j]+M[i,i])*(M[i,j]+M[j,j])) # Covariance between at-site skews. Eq 24 and 25
        }
      }
    }
  }
  
  ## Baseline OLS sigma regression. Eq 15.
  Omega <- diag(nrow(X)) # OLS weighting matrix (identity)
  B.SigReg <- solve(t(X)%*%solve(Omega)%*%X)%*%t(X)%*%solve(Omega)%*%LP3$S # OLS estimated coefficients for k-variable model of LP3 standard deviation.
  Yhat.SigReg <- X%*%B.SigReg # Estimates from sigma regression
  
  ## NOTE: This iteration procedure is currently implemented in WREG v1.05, though not described in the manual.
  ##        It searches across 30 possible values of model-error variance, looks for a sign-change, and then repeats thirty searches on the identified interval.
  ##        It may be possible to improve performance by stopping th loop after a sign change rather than searching the entire space.
  
  ## Set variables to control iterative procedure
  Target<-length(Y) - ncol(X) # Target value of iterations. RHS of Eq 21.
  GInt <- 30 # Number of intervals to consider (doubled below)
  Gstep <- 2*var(Y)/(GInt-1) # Length of step
  ## Coarse intervals.
  Iterations <- matrix(0,ncol=2,nrow=GInt) # Empty dataframe to store results
  Iterations <- data.frame(Iterations); names(Iterations) <- c('GSQ','Deviation') # Formatting empty dataframe.  Will contain estimated model-error variance (GSQ) and deviation from target.
  for (w in 1:GInt) {
    GSQ <- (w-1)*Gstep; Iterations$GSQ[w]<-GSQ; # guess at model error variance
    for (i in 1:length(Y)) {
      for (j in 1:length(Y)) {
        if (i==j) { # Diagonal elements of weighting matrix
          if (SkewAdj) { # if skew-adjusted GLS, use Eq 22 (correct manual to match Giffis and Stedinger (2007) Eq 5; WREG v1.05 code is correct)
            Omega[i,j] <-  GSQ + 
              Yhat.SigReg[i]^2*(1+LP3$K[i]*LP3$G[i]+
                  0.5*LP3$K[i]^2*(1+0.75*LP3$G[i]^2)+
                  Wg[i]*LP3$K[i]*dKdG[i]*(3*LP3$G[i]+0.75*LP3$G[i]^3)+
                  Wg[i]^2*dKdG[i]^2*(6+9*LP3$G[i]^2+1.875*LP3$G[i]^4))/M[i,j] + 
              (1-Wg[i])^2*Yhat.SigReg[i]^2*MSEGR*dKdG[i]^2
          } else { # if normal GLS, use Eq 19.
            Omega[i,j] <-  GSQ + Yhat.SigReg[i]^2*
              (1+LP3$K[i]*LP3$G[i]+0.5*LP3$K[i]^2*
                  (1+0.75*LP3$G[i]^2))/M[i,j]
          }
        } else { # Off-diagonal elements of weighting matrix
          if (SkewAdj) { # if skew-adjusted GLS, use Eq 22 (correct manual to match Giffis and Stedinger (2007) Eq 5; WREG v1.05 code is correct)
            Omega[i,j] <- Rhos[i,j]*Yhat.SigReg[i]*Yhat.SigReg[j]*M[i,j]*
              (1+0.5*LP3$K[i]*LP3$G[i]+0.5*LP3$K[j]*LP3$G[j]+
                  0.5*LP3$K[i]*LP3$K[j]*(Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                  0.5*Wg[j]*LP3$K[i]*LP3$G[j]*dKdG[j]*(3*Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                  0.5*Wg[i]*LP3$K[j]*LP3$G[i]*dKdG[i]*(3*Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                  Wg[i]*Wg[j]*Yhat.SigReg[i]*Yhat.SigReg[j]*dKdG[i]*dKdG[j]*Covgg[i,j])/M[i,i]/M[j,j]
          } else { # if normal GLS, use Eq 19.
            Omega[i,j] <- Rhos[i,j]*Yhat.SigReg[i]*Yhat.SigReg[j]*M[i,j]*
              (1+0.5*LP3$K[i]*LP3$G[i]+0.5*LP3$K[j]*LP3$G[j]+
                  0.5*LP3$K[i]*LP3$K[j]*(Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j]))/M[i,i]/M[j,j]
          }
        }
      }
    }
    B_hat <- solve(t(X)%*%solve(Omega)%*%X)%*%t(X)%*%solve(Omega)%*%Y # Use weighting matrix to estimate regression coefficients
    Iterations$Deviation[w] <- t(Y-X%*%B_hat)%*%solve(Omega)%*%(Y-X%*%B_hat)-Target # Difference between result and target value. Eq 21.
  }
  ## Finer intervals.  Expands iterval with sign change to get closer to zero.
  Signs <- sign(Iterations$Deviation); LastPos <- which(diff(Signs)!=0) # Finds where sign changes from positive to negative
  if (length(LastPos)==0) { # If no change in sign of deviation
    BestPos <- which(abs(Iterations$Deviation)==min(abs(Iterations$Deviation))) # Use the minimum deviation
    GSQ <- Iterations$GSQ[BestPos] # Best estiamte of model-error variance
  } else { # There is a sign change, so expand the interval with sign change.
    Iterations2 <- matrix(0,ncol=2,nrow=GInt) # Empty matrix to store finer intervals.
    Iterations2 <- data.frame(Iterations); names(Iterations) <- c('GSQ','Deviation') # Formatting empty matrix
    Gstep <- (Iterations$GSQ[LastPos+1]-Iterations$GSQ[LastPos])/(GInt-1) # Step for finer intervals
    for (w in 1:GInt) {
      GSQ <- Iterations$GSQ[LastPos]+(w-1)*Gstep; Iterations2$GSQ[w]<-GSQ; # estimate of model-error variance
      for (i in 1:length(Y)) {
        for (j in 1:length(Y)) {
          if (i==j) { # Diagonal elements of weighting matrix
            if (SkewAdj) { # if skew-adjusted GLS, use Eq 22 (correct manual to match Giffis and Stedinger (2007) Eq 5; WREG v1.05 code is correct)
              Omega[i,j] <-  GSQ + 
                Yhat.SigReg[i]^2*(1+LP3$K[i]*LP3$G[i]+
                    0.5*LP3$K[i]^2*(1+0.75*LP3$G[i]^2)+
                    Wg[i]*LP3$K[i]*dKdG[i]*(3*LP3$G[i]+0.75*LP3$G[i]^3)+
                    Wg[i]^2*dKdG[i]^2*(6+9*LP3$G[i]^2+1.875*LP3$G[i]^4))/M[i,j] + 
                (1-Wg[i])^2*Yhat.SigReg[i]^2*MSEGR*dKdG[i]^2
            } else { # if normal GLS, use Eq 19.
              Omega[i,j] <-  GSQ + Yhat.SigReg[i]^2*
                (1+LP3$K[i]*LP3$G[i]+0.5*LP3$K[i]^2*
                    (1+0.75*LP3$G[i]^2))/M[i,j]
            }
          } else { # Off-diagonal elements of weighting matrix
            if (SkewAdj) { # if skew-adjusted GLS, use Eq 22 (correct manual to match Giffis and Stedinger (2007) Eq 5; WREG v1.05 code is correct)
              Omega[i,j] <- Rhos[i,j]*Yhat.SigReg[i]*Yhat.SigReg[j]*M[i,j]*
                (1+0.5*LP3$K[i]*LP3$G[i]+0.5*LP3$K[j]*LP3$G[j]+
                    0.5*LP3$K[i]*LP3$K[j]*(Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                    0.5*Wg[j]*LP3$K[i]*LP3$G[j]*dKdG[j]*(3*Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                    0.5*Wg[i]*LP3$K[j]*LP3$G[i]*dKdG[i]*(3*Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                    Wg[i]*Wg[j]*Yhat.SigReg[i]*Yhat.SigReg[j]*dKdG[i]*dKdG[j]*Covgg[i,j])/M[i,i]/M[j,j]
            } else { # if normal GLS, use Eq 19.
              Omega[i,j] <- Rhos[i,j]*Yhat.SigReg[i]*Yhat.SigReg[j]*M[i,j]*
                (1+0.5*LP3$K[i]*LP3$G[i]+0.5*LP3$K[j]*LP3$G[j]+
                    0.5*LP3$K[i]*LP3$K[j]*(Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j]))/M[i,i]/M[j,j]
            }
          }
        }
      }
      B_hat <- solve(t(X)%*%solve(Omega)%*%X)%*%t(X)%*%solve(Omega)%*%Y # Use weighting matrix to estimate regression coefficients
      Iterations2$Deviation[w] <- t(Y-X%*%B_hat)%*%solve(Omega)%*%(Y-X%*%B_hat)-Target # Difference between result and target value. Eq 21.
    }
    BestPos <- which(Iterations2$Deviation==min(Iterations2$Deviation[Iterations2$Deviation>0])) # Find minimum positive deviation from Eq 21.
    GSQ <- Iterations2$GSQ[BestPos] # best estimate of model-error variance
  }
  ## Calculate Final Omega (weighting matrix)
  for (i in 1:length(Y)) {
    for (j in 1:length(Y)) {
      if (i==j) { # Diagonal elements of weighting matrix
        if (SkewAdj) { # if skew-adjusted GLS, use Eq 22 (correct manual to match Giffis and Stedinger (2007) Eq 5; WREG v1.05 code is correct)
          Omega[i,j] <-  GSQ + 
            Yhat.SigReg[i]^2*(1+LP3$K[i]*LP3$G[i]+
                0.5*LP3$K[i]^2*(1+0.75*LP3$G[i]^2)+
                Wg[i]*LP3$K[i]*dKdG[i]*(3*LP3$G[i]+0.75*LP3$G[i]^3)+
                Wg[i]^2*dKdG[i]^2*(6+9*LP3$G[i]^2+1.875*LP3$G[i]^4))/M[i,j] + 
            (1-Wg[i])^2*Yhat.SigReg[i]^2*MSEGR*dKdG[i]^2
        } else { # if normal GLS, use Eq 19.
          Omega[i,j] <-  GSQ + Yhat.SigReg[i]^2*
            (1+LP3$K[i]*LP3$G[i]+0.5*LP3$K[i]^2*
                (1+0.75*LP3$G[i]^2))/M[i,j]
        }
      } else { # Off-diagonal elements of weighting matrix
        if (SkewAdj) { # if skew-adjusted GLS, use Eq 22 (correct manual to match Giffis and Stedinger (2007) Eq 5; WREG v1.05 code is correct)
          Omega[i,j] <- Rhos[i,j]*Yhat.SigReg[i]*Yhat.SigReg[j]*M[i,j]*
            (1+0.5*LP3$K[i]*LP3$G[i]+0.5*LP3$K[j]*LP3$G[j]+
                0.5*LP3$K[i]*LP3$K[j]*(Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                0.5*Wg[j]*LP3$K[i]*LP3$G[j]*dKdG[j]*(3*Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                0.5*Wg[i]*LP3$K[j]*LP3$G[i]*dKdG[i]*(3*Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                Wg[i]*Wg[j]*Yhat.SigReg[i]*Yhat.SigReg[j]*dKdG[i]*dKdG[j]*Covgg[i,j])/M[i,i]/M[j,j]
        } else { # if normal GLS, use Eq 19.
          Omega[i,j] <- Rhos[i,j]*Yhat.SigReg[i]*Yhat.SigReg[j]*M[i,j]*
            (1+0.5*LP3$K[i]*LP3$G[i]+0.5*LP3$K[j]*LP3$G[j]+
                0.5*LP3$K[i]*LP3$K[j]*(Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j]))/M[i,i]/M[j,j]
        }
      }
    }
  }
  ## Control output
  GLS.Weights <- list(GSQ=GSQ,Omega=Omega) # Output contains model-error variance and weighting matrix.
  return(GLS.Weights)
}

#' Distance calculation (WREG)
#'
#'@description
#' \code{Dist.WREG} calculates the distance between two points defined
#'  by Latitude-Longitude coordinates in decimal degrees.
#'
#'@param Lat1 Latitude of the first point, in decimal degrees.
#'@param Long1 Lonigtude of the first point, in decimal degrees.
#'@param Lat2 Latitude of the second point, in decimal degrees.
#'@param Long2 Longitude of the second point, in deceimal degrees.
#'@param method Idicates which technique to use for distance calculation.  See details.
#'
#'@details
#'\code{Dist.WREG} is capable of using two techniques to calculate intersite distances.  
#'\code{method==1} indicates that the "Nautical Mile" approximation should be used.
#'This is the function that is currently employed by WREG v. 1.05.  Each arcminute is equal to 1852 meters.
#'\code{method==2} indicates that the Haversine approximation should be used.
#'
#'@return Returns the distance between the two sites in miles.
#'@export
Dist.WREG <- function(Lat1,Long1,Lat2,Long2,method=c(1,2)) {
  # William Farmer, USGS, January 23, 2015
  
  if (method==1) {
    ## Nautical mile conversion
    Dist <- sqrt((abs(Lat2-Lat1)*1852*60)^2+(abs(Long2-Long1)*1852*60)^2)*0.6214/1000 # Intersite distance, miles
  } else if (method==2) {
    R <- 6371/1.609 # Radius of the earth in miles
    ## Convert Lat/Long to Radians
    Lat1 <- Lat1*pi/180
    Long1 <- Long1*pi/180
    Lat2 <- Lat2*pi/180
    Long2 <- Long2*pi/180
    ## Haversine formula
    a <- (sin(0.5*(Lat2-Lat1)))^2+cos(Lat1)*cos(Lat2)*(sin(0.5*(Long2-Long1))^2)
    c <- 2*atan2(sqrt(a),sqrt(1-a))
    Dist <- R*c # Intersite distance, miles
  }
  return(Dist)
}

#' Region-of-Influence Regression (WREG)
#'
#'@description
#'\code{WREG.ROI} implements region-of-influence regression in the WREG framework.
#'
#' @param Y The dependent variable of interest, with any transformations already 
#' applied.
#' @param X The independent variables in the regression, with any transformations 
#' already applied.  Each row represents a site and each column represents
#' a particular independe variable.  (If a leading constant is used, it should be 
#' included here as a leading column of ones.)  The rows must be in the same order as
#' the dependent variables in \code{Y}.
#' @param Reg A string indicating which type of regression should be applied.
#' The options include: \dQuote{OLS} for ordinary least-squares regression,
#' \dQuote{WLS} for weighted least-squares regression, \dQuote{GLS} for generalized
#' least-squares regression, with no uncertainty in regional skew and \dQuote{GLSskew}
#' for generalized least-squares regression with uncertainty in regional skew.  
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
#' @param ROI A string indicating how to define the region of influence.  \dQuote{PRoI}
#' signifies physiographic, independent or predictor-variable region of influence.
#' \dQuote{GRoI} calls for a geographic region of influence. \dQuote{HRoI} requests
#' a hybrid region of influence.  Details on each approach are provided in the manual
#' for WREG v. 1.0.
#' @param n The number of sites to include in the region of influence.
#' @param D Required for \dQuote{HRoI}, the geographic limit within which to search for a 
#' physiographic region of influence.  In WREG v. 1.05 (see \code{Legacy} below),
#' this is interpretted in meters.  Elsewise, this is interpretted as miles.
#' @param DistMeth Required for \dQuote{GLS} and \dQuote{GLSskew}.  A value of \code{1} 
#' indicates that the "Nautical Mile" approximation should be used to calculate inter-site
#' distances.  A value of \code{2} designates the Haversine approximation.  See 
#' \code{\link{Dist.WREG}}.  The default value is \code{2}.  (See the \code{Legacy} 
#' details below.)
#' @param Legacy A logical.  A value of \code{TRUE} forces the WREG program to behave 
#' identically to WREG v. 1.05, with BUGS and all.  It will force \code{TY=2} and
#' \code{DistMeth=1}.  For ROI regressions, it will also require a specific calculation 
#' for weighing matrices in \dQuote{WLS} (\code{\link{Omega.WLS.ROImatchMatLab}}), 
#' \dQuote{GLS}, and \dQuote{GLSskew} (see \code{\link{Omega.GLS.ROImatchMatLab}}).
#' \code{Legacy} also forces the distance limit \code{D} to be interpretted in meters.
#'
#'@details The support for region-of-influence regression is described in the manual
#'of WREG v. 1.0.  \code{WREG.RoI} iterates through the sites of \code{Y}, defines a 
#'region of influence and implements the specified regression by calling 
#'\code{\link{WREG.MLR}}.
#'
#' The logical handle \code{Legacy} has been included to test that this program returns the 
#' same results as WREG v. 1.05.  In the development of this code, some idiosyncrasies of 
#' the MatLab code for WREG v. 1.05 became apparent.  Setting \code{Legacy} equal to 
#' \code{TRUE} forces the code to use the same idiosycrasies as WREG v. 1.05.  Some of 
#' these idosyncrasies may be bugs in the code.  Further analysis is needed.  For 
#' information on the specific idiosyncrasies, see the notes for the \code{Legacy} 
#' input and the links to other functions in this package.
#'
#'@return As with \code{\link{WREG.MLR}}, \code{WREG.RoI} returns a large list of 
#'regression parameters and metrics.  This list varies depending on the \code{Reg}
#'specified, but may contain:
#' \item{fitted.values}{A vector of model estimates from the regression model.}
#' \item{residuals}{A vector of model residuals.}
#' \item{PerformanceMetrics}{A list of four elements.  These represent approximate
#' performance regression across all of the region-of-influence regressions. These 
#' include the mean squared error of residuals (\code{MSE}), the coefficient of
#' determination (\code{R2}), the adjusted coefficient of determination (\code{R2_adj})
#' and the root mean squared error (\code{RMSE}, in percent).  Details on the 
#' appropriateness and applicability of performance metrics can be found in the WREG manual.}
#' \item{Coefficients}{A list composed of four elements: (1) \code{Values} contains 
#' the regression coefficeints estimated for the model built around each observation, 
#' (2) \code{\sQuote{StanError}} contains the standard errors of each regression coefficient 
#' for the ROI regressions around each observations, (3) \code{TStatistic} 
#' contains the Student's T-statistic of each regression coefficient for the ROI 
#' regression built around each observation and (4) \code{pValue} contains the significance 
#' probability of each regression coefficient for the ROI regressions built around each
#' observation.  Each element of the list is a matrix the same size as \code{X}}
#' \item{ROI.Regressions}{A list of elements and outputs from each individual ROI
#' regression.  These include a matrix of the sites used in each ROI regression 
#' (\code{Sites.Used}, a \code{length(Y)}-by-\code{n} matrix), a matrix of the 
#' geographic distances between the selected sites in each ROI regression (\code{
#' Gdist.Used}, a \code{length(Y)}-by-\code{n} matrix), a matrix of the 
#' physiographic distances between the selected sites in each ROI regression (\code{
#' Pdist.Used}, a \code{length(Y)}-by-\code{n} matrix), a matrix of the 
#' observations used in each ROI regression (\code{Obs.Used}, a \code{length(Y)}-by-
#' \code{n} matrix), a matrix of model fits in the region of influence (\code{Fits}, 
#' a \code{length(Y)}-by-\code{n} matrix), a matrix of model residuals in the region of 
#' influence (\code{Residuals}, a \code{length(Y)}-by-\code{n} matrix), a matrix of 
#' leverages for each ROI regression (\code{Leverage}, a \code{length(Y)}-by-\code{n} 
#' matrix), a matrix of influences for each ROI regression (\code{Influence}, a \code{
#' length(Y)}-by-\code{n} matrix), a logical matrix indicating if the leverage is 
#' significant in each ROI regression (\code{Leverage.Significance}, a \code{length(Y)}
#' -by-\code{n} matrix), a logical matrix indicating if the influence is 
#' significant in each ROI regression (\code{Influence.Significance}, a \code{length(Y)}
#' -by-\code{n} matrix), a vector of critical leverage values for each ROI regression
#' (\code{Leverage.Limits}), a vector of critical influence values for each ROI
#' regression (\code{Influence.Limits}) and a list of performance metrics for each 
#' ROI regression (\code{PerformanceMetrics}).  The last element, 
#' \code{PerformanceMetrics} is identical to the same output from \code{\link{
#' WREG.MLR}} excpet that every element is multiplied by the number of observations
#' so as to capture the individual performance of each ROI regression.}
#' \item{ROI.InputParameters}{A list of input parameters to record the controls on 
#' the ROI regression.  \code{D} idicates the limit used in \dQuote{HRoI}.
#' \code{n} indicates the size of the region of influence.  \code{ROI} is a string 
#' indicating the type of region of influence.  \code{Legacy} is a logical indicating
#' if the WREG v. 1.05 idiosycrasies were implemented.}
#'@export
WREG.RoI <- function(Y,X,Reg=c('OLS','WLS','GLS','GLSskew'),RecordLengths=NA,LP3=NA,
  alpha=0.01,theta=0.98,BasinChars=NA,MSEGR=NA,TY=2,Peak=T,
  ROI=c('PRoI','GRoI','HRoI'),n,D=250,DistMeth=2,Legacy=FALSE) {
  # William Farmer, USGS, January 23, 2015
  
  ## Add controls to meet Legacy demands
  ##    NOTE: Legacy forces program to return the same results as WREG v 1.05.
  if (Legacy) { # If legacy is indicated, override custom inputs.
    TY <- 2 # WREG v1.05 does not read this input correctly.
    DistMeth <- 1 # WREG v1.05 uses "Nautical Mile" approximation
  }
  
  ## Empty vectors for output
  Sites.Used <- matrix(NA,nrow=length(Y),ncol=n) # Empty matrix.  Sites used for region of influence.
  Gdist.Used <- matrix(NA,nrow=length(Y),ncol=n) # Empty matrix.  Geographic distance to sites used for region of influence.
  Pdist.Used <- matrix(NA,nrow=length(Y),ncol=n) # Empty matrix.  Physiographic distance to sites used for region of influence.
  Obs.Used <- matrix(NA,nrow=length(Y),ncol=n) # Empty matrix.  Observed predictands at sites used for region of influence.
  Fits <- matrix(NA,nrow=length(Y),ncol=n) # Empty matrix.  Model fits at sites used for region of influence.
  Residuals <- matrix(NA,nrow=length(Y),ncol=n) # Empty matrix.  Model residuals at sites used for region of influence.
  Leverage.v <- matrix(NA,nrow=length(Y),ncol=n) # Empty matrix.  Leverage of sites used for region of influence.
  Influence.v <- matrix(NA,nrow=length(Y),ncol=n) # Empty matrix.  Influence of sites used for region of influence.
  Leverage.Significance <- matrix(NA,nrow=length(Y),ncol=n) # Empty matrix.  Logical.  Sites used for region of influence with significant leverage.
  Influence.Significance <- matrix(NA,nrow=length(Y),ncol=n) # Empty matrix.  Logical.  Sites used for region of influence with significant influence.
  Leverage.Limit <- vector(length=length(Y)) # Empty vector.  Critical values of leverage.
  Influence.Limit <- vector(length=length(Y)) # Empty vector.  Critical values of influence.
  Coef.Values <- matrix(NA,nrow=length(Y),ncol=ncol(X)) # Empty matrix.  Coefficents estimated from region of influence.
  Coef.SE <- matrix(NA,nrow=length(Y),ncol=ncol(X)) # Empty matrix.  Standard errors of coefficients estimated from region of influence.
  Coef.T <- matrix(NA,nrow=length(Y),ncol=ncol(X)) # Empty matrix.  T-statistics of coefficients estimated from region of influence.
  Coef.p <- matrix(NA,nrow=length(Y),ncol=ncol(X)) # Empty matrix.  p-values of coefficients estimated from region of influence.
  Y_hat <- vector(length=length(Y)) # Empty vector.  Region-of-Influence model fits.
  ROI.Omegas <- array(dim=c(length(Y),n,n)) # Empty array.  Weighting matrices from region-of-influence regressions.
  ROI.PerfMets <- list(MSE=vector(length=length(Y)),R2=vector(length=length(Y)),R2_adj=vector(length=length(Y)),RMSE=vector(length=length(Y))) # Performance metrics for output (basic, for OLS)
  if (is.element(Reg,c('WLS','GLS','GLSskew'))) { # non-OLS requires additional performance metrics
    ROI.PerfMets$R2_pseudo <- vector(length=length(Y)) # Empty vector.  Pseudo coefficient of determination of the individual ROI regression
    ROI.PerfMets$AVP=vector(length=length(Y)) # Empty vector.  Average variance of prediction for the individual ROI regression
    ROI.PerfMets$Sp=vector(length=length(Y)) # Empty vector.  Standard error of prediction, in percent, for the individual ROI regression
    ROI.PerfMets$VP.PredVar=matrix(NA,nrow=length(Y),ncol=n) # Empty vector.  Individual variances of prediction for the individual ROI regression
    ROI.PerfMets$ModErrVar=vector(length=length(Y)) # Empty vector.  Model-error variance of individual ROI regression
    ROI.PerfMets$StanModErr=vector(length=length(Y)) # Empty vector.  Standardized model-error variance, in percent, of individual ROI regression
  }
  
  ## Determine standard deviations of explanatory variables
  if (length(unique(X[,1]))==1) { # If the first predictor is a constant
    X.nocon <- X[,2:ncol(X)] # remove leading constant from model
  } else {
    X.nocon <- X # there is no leading constant
  }
  SDs <- t(matrix(rep(apply(X.nocon,2,sd),times=length(Y)),ncol=length(Y))) # Standard deviations of each predictor.
  
  ## Cycle through sites.  For each site, determine region, perform regression, save output.
  for (i in 1:length(Y)) { # Loop through sites...
    x0.i <- X[i,] # Predictors at target site
    x0.i.nocon <- t(matrix(rep(X.nocon[i,],times=length(Y)),ncol=length(Y))) # Predictors at target site, less a constant
    ### Physiographic Euclidian distance
    Pdist <- sqrt(rowSums(((x0.i.nocon-X.nocon)/SDs)^2)) # Euclidian distance in independent-variable space.  Eq 47.
    Pdist[i] <- Inf # To block self identification.
    ### Geographic distance distance
    Gdist <- vector(length=length(Y)) # Empty vector for geographic distances
    for (j in 1:length(Y)) {
      if (i!=j) {
        #### Geographic distance
        Gdist[j] <- Dist.WREG(Lat1=BasinChars$Lat[i],Long1=BasinChars$Long[i],Lat2=BasinChars$Lat[j],Long2=BasinChars$Long[j],method=DistMeth) # Intersite distance, miles
      } else {
        Gdist[j] <- Inf # To block self identification.
      }
    }
    if (Legacy) {
      Gdist <- Gdist*1000/0.6214 # Original matlab code does not specify the units for the limit.  The code implies meters.  Dist.WREG returns miles.
      # BUG: Should correct this so that D is interpretted in miles.
    }
    
    ### Identify region of influence
    if (ROI=='PRoI') {
      #### Physiographic Region-of-Influence
      temp <- sort.int(Pdist,index.return=T)
      NDX <- temp$ix[1:n] # Sites to use in this regression
    } else if (ROI=='GRoI') {
      #### Geographic Region-of-Influence
      temp <- sort.int(Gdist,index.return=T)
      NDX <- temp$ix[1:n] # Sites to use in this regression
    } else if (ROI=='HRoI') {
      #### Hybrid Region-of-Influence
      if (sum(Gdist<=D)<n) { # Not enough sites for hybrid selection; default to GRoI
        temp <- sort.int(Gdist,index.return=T)
        NDX <- temp$ix[1:n] # Sites to use in this regression
      } else { # Conduct hybrid ranking
        ProximateSites <- which(Gdist<=D,arr.ind=T)
        temp <- sort.int(Pdist[ProximateSites],index.return=T)
        NDX <- ProximateSites[temp$ix[1:n]] # Sites to use in this regression
      }
    }
    
    ### Perform regressions
    Y.i <- Y[NDX] # Predictands from region of influence
    X.i <- X[NDX,] # Predictors from region of influence
    if (!is.na(sum(c(RecordLengths)))&&is.matrix(RecordLengths)) {# GLS
      RecordLengths.i <- RecordLengths[NDX,NDX] # Record lengths from region of influence
    } else if (!is.na(sum(c(RecordLengths)))&&is.vector(RecordLengths)) { # WLS
      RecordLengths.i <- RecordLengths[NDX] # Record lengths from region of influence
    }
    BasinChars.i <- BasinChars[NDX,] # Basin characteristics (IDs, Lat, Long) from region of influence.
    LP3.i <- data.frame(LP3)[NDX,] # LP3 parameters from region of influence
    if (Legacy) { # Apply regressions to match MatLab WREG v 1.05.
      if (Reg=='WLS') { # Use subroutine to correct WLS weights to meet v1.05 idiosyncrasies
        WeightFix <- Omega.WLS.ROImatchMatLab(Y.all=Y,X.all=X,LP3.all=LP3,RecordLengths.all=RecordLengths,NDX=NDX)
        Reg.i <- WREG.MLR(Y=Y.i,X=X.i,x0=x0.i,Reg='CustomWeight',RecordLengths=RecordLengths.i,LP3=LP3.i,CustomWeight=WeightFix,Legacy=Legacy)
      } else if (is.element(Reg,c('GLS','GLSskew'))) { # Use subroutine to correct GLS and GLSskew weights to meet v1.05 idiosyncrasies
        # "Corrected" weighting matrix to match MATLAB code.  Returns weighting matrix and var.moderror.k
        WeightFix.k <- Omega.GLS.ROImatchMatLab(alpha=alpha,theta=theta,Independent=BasinChars.i,X=X.i,Y=Y.i,RecordLengths=RecordLengths.i,LP3=LP3.i,MSEGR=MSEGR,TY=TY,Peak=Peak,X.all=X,LP3.all=LP3,DistMeth=DistMeth)
        # "Corrected" to match MATLAB code.  Returns weighting matrix and var.moderror.0
        WeightFix.0 <- Omega.GLS.ROImatchMatLab(alpha=alpha,theta=theta,Independent=BasinChars.i,X=matrix(1,ncol=1,nrow=nrow(X.i)),Y=Y.i,RecordLengths=RecordLengths.i,LP3=LP3.i,MSEGR=MSEGR,TY=TY,Peak=Peak,X.all=matrix(1,ncol=1,nrow=nrow(X)),LP3.all=LP3,DistMeth=DistMeth)
        WeightFix <- list(Omega=WeightFix.k$Omega,var.modelerror.0=WeightFix.0$GSQ,var.modelerror.k=WeightFix.k$GSQ)
        Reg.i <- WREG.MLR(Y=Y.i,X=X.i,x0=x0.i,Reg='CustomWeight',CustomWeight=WeightFix,Legacy=Legacy) 
      } else if (Reg=='OLS') { # Apply OLS normally.
        Reg.i <- WREG.MLR(Y=Y.i,X=X.i,x0=x0.i,Reg=Reg,Legacy=Legacy)
      }
    } else { # Aply WREG with "bugs" corrected.
      Reg.i <- WREG.MLR(Y=Y.i,X=X.i,x0=x0.i,Reg=Reg,RecordLengths=RecordLengths.i,LP3=LP3.i,alpha=alpha,theta=theta,BasinChars=BasinChars.i,MSEGR=MSEGR,TY=TY,Peak=Peak,DistMeth=DistMeth,Legacy=Legacy)
    }
    
    ### Store outputs
    Y_hat[i] <- Reg.i$Y.ROI # Region-of-Influence model fits.
    Sites.Used[i,] <- NDX # Sites used for region of influence.
    Gdist.Used[i,] <- Gdist[NDX] # Geographic distance to sites used for region of influence.
    Pdist.Used[i,] <- Pdist[NDX] # Physiographic distance to sites used for region of influence.
    Obs.Used[i,] <- Y.i # Observed predictands at sites used for region of influence.
    Fits[i,] <- Reg.i$fitted.values # Model fits at sites used for region of influence.
    Residuals[i,] <- Reg.i$residuals # Model residuals at sites used for region of influence.
    Leverage.v[i,] <- Reg.i$ResLevInf$Leverage # Leverage of sites used for region of influence.
    Influence.v[i,] <- Reg.i$ResLevInf$Influence # Influence of sites used for region of influence.
    Leverage.Significance[i,] <- Reg.i$LevInf.Sig[,1] # Logical.  Sites used for region of influence with significant leverage.
    Influence.Significance[i,] <- Reg.i$LevInf.Sig[,2] # Logical.  Sites used for region of influence with significant influence.
    Leverage.Limit[i] <- Reg.i$LevLim  # Critical values of leverage.
    Influence.Limit[i] <- Reg.i$InflLim  # Critical values of influence.
    Coef.Values[i,] <- Reg.i$Coefs$Coefficient  # Coefficents estimated from region of influence.
    Coef.SE[i,] <- Reg.i$Coefs$'Standard Error' # Standard errors of coefficients estimated from region of influence.
    Coef.T[i,] <- Reg.i$Coefs$tStatistic # T-statistics of coefficients estimated from region of influence.
    Coef.p[i,] <- Reg.i$Coefs$pValue # p-values of coefficients estimated from region of influence.
    ROI.Omegas[i,,] <- Reg.i$Weighting # Weighting matrices from region-of-influence regressions.
    ROI.PerfMets$MSE[i] <- Reg.i$PerformanceMetrics$MSE # Mean squared-error of the individual ROI regression
    ROI.PerfMets$R2[i] <- Reg.i$PerformanceMetrics$R2 # Coefficient of determination of the individual ROI regression
    ROI.PerfMets$R2_adj[i] <- Reg.i$PerformanceMetrics$R2_adj # Adjusted coefficient of determination of the individual ROI regression
    ROI.PerfMets$RMSE[i] <- Reg.i$PerformanceMetrics$RMSE # Root-mean-squared error of the individual ROI regression
    if (is.element(Reg,c('WLS','GLS','GLSskew'))) { # non-OLS requires additional performance metrics
      ROI.PerfMets$R2_pseudo[i] <- Reg.i$PerformanceMetrics$R2_pseudo # Pseudo coefficient of determination of the individual ROI regression
      ROI.PerfMets$AVP[i] <- Reg.i$PerformanceMetrics$AVP # Average variance of prediction for the individual ROI regression
      ROI.PerfMets$Sp[i] <- Reg.i$PerformanceMetrics$Sp # Standard error of prediction, in percent, for the individual ROI regression
      ROI.PerfMets$VP.PredVar[i,] <- Reg.i$PerformanceMetrics$VP.PredVar # Individual variances of prediction for the individual ROI regression
      ROI.PerfMets$ModErrVar[i] <- Reg.i$PerformanceMetrics$ModErrVar # Model-error variance of individual ROI regression
      ROI.PerfMets$StanModErr[i] <- Reg.i$PerformanceMetrics$StanModErr # Standardized model-error variance, in percent, of individual ROI regression
    }
  }
  
  e <- Y-Y_hat # ROI model residuals
  MSE <- sum(e^2)/(nrow(X)-ncol(X)) # Mean square-error of ROI estimates
  if (Legacy) MSE <- mean(e^2) # BUG: ROI code takes mean of squared residuals without accounting for parameters.
  SSR <- sum(e^2) # Residual sum of squares from ROI estimates
  SST <- sum((Y-mean(Y))^2) # Total sum of squares from observations
  R2 <- 1 - SSR/SST # Coefficient of determination of overall ROI estimates
  R2_adj <- 1 - SSR*(nrow(X)-1)/SST/(nrow(X)-ncol(X)) # Adjusted coefficient of determination of overall ROI estimates
  RMSE <-100*sqrt(exp(log(10)*log(10)*MSE)-1) # Root-mean-squared error, in percent, of overall ROI estimates
  PerfMet <- list(MSE=MSE,R2=R2,R2_adj=R2_adj,RMSE=RMSE) # Performance metrics (for output)
  
  
  Output.ROI <- list(fitted.values=Y_hat,residuals=e,
    PerformanceMetrics=PerfMet,
    Coefficients=list(Values=Coef.Values,StanError=Coef.SE,TStatistic=Coef.T,pValue=Coef.p),
    ROI.Regressions=list(Sites.Used=Sites.Used,
      Gdist.Used=Gdist.Used,Pdist.Used=Pdist.Used,
      Obs.Used=Obs.Used,Fits=Fits,Residuals=Residuals,
      Leverage=Leverage.v,Influence=Influence.v,
      Leverage.Significance=Leverage.Significance,
      Influence.Significance=Influence.Significance,
      Leverage.Limits=Leverage.Limit,Influence.Limits=Influence.Limit,
      PerformanceMetrics=ROI.PerfMets),
    ROI.InputParameters=list(D=D,n=n,ROI=ROI,Legacy=Legacy))
  return(Output.ROI)
}

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

#' Weighing Matrix for ROI-GLS/skew (WREG)
#'
#'@description
#'\code{Omega.GLS.ROImatchMatLab} calculates the weighting function for a 
#'generalized least-squares regression using regions-of-influence.  This is largely
#'legacy code to match WREG v. 1.05 idiosyncrasies.
#'
#' @param alpha A number, required only for \dQuote{GLS} and \dQuote{GLSskew}.  
#' \code{alpha} is a parameter used in the estimated cross-correlation between site
#' records.  See equation 20 in the WREG v. 1.05 manual.  The arbitrary, default value 
#' is 0.01.  The user should fit a different value as needed.
#' @param theta A number, required only for \dQuote{GLS} and \dQuote{GLSskew}.  
#' \code{theta} is a parameter used in the estimated cross-correlation between site
#' records.  See equation 20 in the WREG v. 1.05 manual.  The arbitrary, default value 
#' is 0.98.  The user should fit a different value as needed.
#' @param Independent A dataframe containing three variables: \code{StationID} is the 
#' numerical identifier (without a leading zero) of each site, \code{Lat} is the latitude
#' of the site, in decimal degrees, and \code{Long} is the longitude of the site, in decimal
#' degrees.  The sites must be presented in the same order as \code{Y}.  Required only for
#' \dQuote{GLS} and \dQuote{GLSskew}.  This includes only sites in the current region 
#' of influence.
#' @param X The independent variables in the regression, with any transformations 
#' already applied.  Each row represents a site and each column represents
#' a particular independe variable.  (If a leading constant is used, it should be 
#' included here as a leading column of ones.)  The rows must be in the same order as
#' the dependent variables in \code{Y}.  This includes only sites in the current region 
#' of influence.
#' @param Y The dependent variable of interest, with any transformations already 
#' applied.  This includes only sites in the current region of influence.
#' @param RecordLengths This input is required for \dQuote{WLS}, \dQuote{GLS} and 
#' \dQuote{GLSskew}.  For \dQuote{GLS} and \dQuote{GLSskew}, \code{RecordLengths} 
#' should be a matrix whose rows and columns are in the same order as \code{Y}.  Each 
#' \code{(r,c)} element represents the length of concurrent record between sites 
#' \code{r} and \code{c}.  The diagonal elements therefore represent each site's full
#' record length.  For \dQuote{WLS}, the only the at-site record lengths are needed.
#' In the case of \dQuote{WLS}, \code{RecordLengths} can be a vector or the matrix 
#' described for \dQuote{GLS} and \dQuote{GLSskew}.  This includes only sites in the 
#' current region of influence.
#' @param LP3 A dataframe containing the fitted Log-Pearson Type III standard
#' deviate, standard deviation and skew for each site.  The names of this data frame are
#' \code{S}, \code{K} and \code{G}.  For \dQuote{GLSskew}, the regional skew value must 
#' also be provided in a variable called \code{GR}.  The order of the rows must be the same
#' as \code{Y}.  This includes only sites in the current region of influence.
#' @param MSEGR A number. The mean squared error of the regional skew.  Required only for
#' \dQuote{GLSskew}.
#' @param TY A number.  The return period of the event being modeled.  Required only for 
#' \dQuote{GLSskew}.  The default value is \code{2}.  (See the \code{Legacy} details below.)
#' @param Peak A logical.  Indicates if the event being modeled is a peak flow event
#' or a low-flow event.  \code{TRUE} indicates a peak flow, while \code{FALSE} indicates
#' a low-flow event.
#' @param X.all The independent variables for all sites in the network, with any 
#' transformations already applied.  Each row represents a site and each column 
#' represents a particular independe variable.  (If a leading constant is used, it 
#' should be included here as a leading column of ones.)
#' @param LP3.all A dataframe containing the fitted Log-Pearson Type III standard
#' deviate, standard deviation and skew for all sites in the network.  The names of 
#' this data frame are \code{S}, \code{K} and \code{G}.  For \dQuote{GLSskew}, the 
#' regional skew value must also be provided in a variable called \code{GR}.  The 
#' order of the rows must be the same as \code{Y}.
#' @param DistMeth Required for \dQuote{GLS} and \dQuote{GLSskew}.  A value of \code{1} 
#' indicates that the "Nautical Mile" approximation should be used to calculate inter-site
#' distances.  A value of \code{2} designates the Haversine approximation.  See 
#' \code{\link{Dist.WREG}}.  The default value is \code{2}.  (See the \code{Legacy} 
#' details below.)
#'
#'@details This is a legacy function that matches the idiosyncrasies of WREG v. 1.05.
#'This includes using all sites to implement the sigma regression.
#'
#'See \code{\link{Omega.GLS}} for more information on the \dQuote{GLS} 
#'weighting estimates.
#'
#'This function will become obsolete once all idiosyncrasies are assessed.
#'
#'@return \code{Omega.GLS.ROImatchMatLab} returns a list with two elements:
#'\item{GSQ}{The estimated model error variance.}
#'\item{Omega}{The estimated weighting matrix.}
#'@export
Omega.GLS.ROImatchMatLab <- function(alpha=0.01,theta=0.98,Independent,X,Y,RecordLengths,
  LP3,MSEGR=NA,TY=2,Peak=T,X.all,LP3.all,DistMeth=2) {
  # William Farmer, January 27, 2015
  #
  # WREG v 1.05 (MatLab) is inconsistent in which components are used to develop equations in RoI-WLS
  # Eq. 15 is implemented with all sites
  
  ## Determining if skew adjustment is requested
  SkewAdj<-F # default: no skew adjustment
  if (!is.na(MSEGR)) {
    SkewAdj<-T # If user provides a mean squared-error of regional skew, then use skew adjustment.
  }
  
  ## Create distance matrix and concurrent record lengths
  ##    (For skew-adjusted GLS: Also calculates the LP3 partial derivatives, mean squared-errors of at-site skew and the variance of at-site skew.)
  Dists <- matrix(NA,ncol=length(Y),nrow=length(Y)) # Empty matrix for intersite distances
  M <- RecordLengths # Just renaming input for ease.  (Should probably correct later...)
  if (SkewAdj) { # Make empty vectors for GLS-skew
    dKdG <- vector(length=length(Y)) # Empty vector for LP3 partial derivatives
    MSEg <- vector(length=length(Y)) # Empty vector for mean squared-error of at-site skew
    Varg <- vector(length=length(Y)) # Empty vector for variance of at-site skew
    ### Convert return period into probability
    if (Peak) { # if a peak flow is being estimated
      Zp <- -qnorm(1/TY) 
    } else { # if a low flow is being estimated
      Zp <- qnorm(1/TY)
    }
  }
  for (i in 1:length(Y)) {
    for (j in 1:length(Y)) {
      if (i!=j) {
        ### Calculate intersite distance via subroutine.
        ###     DistMeth==1 applies the 'Nautical Mile' approximation from WREG v 1.05
        ###     DistMeth==2 applies Haversine approximation
        Dists[i,j] <- Dist.WREG(Lat1=Independent$Lat[i],Long1=Independent$Long[i],Lat2=Independent$Lat[j],Long2=Independent$Long[j],method=DistMeth) # Intersite distance, miles
      }
    }
    if (SkewAdj) { # Additional calculations for skew-adjusted GLS
      ### LP3 Partial Derivative
      dKdG[i] <- (Zp^2-1)/6+LP3$G[i]*(Zp^3-6*Zp)/54-LP3$G[i]^2*(Zp^2-1)/72+Zp*LP3$G[i]^3/324+5*LP3$G[i]^4/23328 # LP3 partial derivative. Eq 23.
      ### Variance of at-site skew
      a <- -17.75/M[i,i]^2+50.06/M[i,i]^3 # Coefficeint for calculation of the variance of at-site skew.  Eq 27.
      b1 <- 3.92/M[i,i]^0.3-31.1/M[i,i]^0.6+34.86/M[i,i]^0.9  # Coefficeint for calculation of the variance of at-site skew.  Eq 28.
      c1 <- -7.31/M[i,i]^0.59+45.9/M[i,i]^1.18-86.5/M[i,i]^1.77  # Coefficeint for calculation of the variance of at-site skew.  Eq 29.
      Varg[i] <- (6/M[i,i]+a)*(1+LP3$GR[i]^2*(9/6+b1)+LP3$GR[i]^4*(15/48+c1)) # Variance of at-site skew. Eq 26.
      ### Mean squared-error of at-site skew
      ###   These equations are not included or described in v1.05 manual.  They are similar to Eq 28 and 29.
      ###   These equations, which appear in the v1.05 code, are documented in Griffis and Stedinger (2009); Eq 3, 4, 5 and 6.
      b <- 3.93/M[i,i]^0.3-30.97/M[i,i]^0.6+37.1/M[i,i]^0.9 # Coefficient for the calculation of mean squared-error of at-site skew.
      c <- -6.16/M[i,i]^0.56+36.83/M[i,i]^1.12-66.9/M[i,i]^1.68 # Coefficient for the calculation of mean squared-error of at-site skew.
      MSEg[i] <- (6/M[i,i]+a)*(1+LP3$G[i]^2*(9/6+b)+LP3$G[i]^4*(15/48+c)) # Mean squared-error of at-site skew.
    }
  }
  Rhos <- theta^(Dists/(alpha*Dists+1)) # Estimated intersite correlation. Eq 20.
  
  if (SkewAdj) { # if skew-adjusted GLS
    ### Calculate skew weights.
    Wg <- MSEGR/(MSEg+MSEGR) # skew weight. Eq 17.
    ### Calculate covariances between at-site skews
    Covgg <- matrix(NA,ncol=length(Y),nrow=length(Y)) # Empty matrix for covariances between at-site skews.
    for (i in 1:length(Y)) {
      for (j in 1:length(Y)) {
        if (i!=j) {
          Covgg[i,j] <- M[i,j]*sign(Rhos[i,j])*abs(Rhos[i,j])^3*sqrt(Varg[i]*Varg[j])/sqrt((M[i,j]+M[i,i])*(M[i,j]+M[j,j])) # Covariance between at-site skews. Eq 24 and 25
        }
      }
    }
  }
  
  ## Baseline OLS sigma regression. Eq 15.
  Omega <- diag(nrow(X.all)) # OLS weighting matrix (identity)
  B.SigReg <- solve(t(X.all)%*%solve(Omega)%*%X.all)%*%t(X.all)%*%solve(Omega)%*%LP3.all$S # OLS estimated coefficients for k-variable model of LP3 standard deviation.
  # BUG: MatLab fits beta regression (Eq 15) using all sites.  See lines 1305-1316 and 1372-1378 of WREG v 1.05
  Yhat.SigReg <- X%*%B.SigReg # Estimates from sigma regression
  
  ## NOTE: This iteration procedure is currently implemented in WREG v1.05, though not described in the manual.
  ##        It searches across 30 possible values of model-error variance, looks for a sign-change, and then repeats thirty searches on the identified interval.
  ##        It may be possible to improve performance by stopping th loop after a sign change rather than searching the entire space.
  
  ## Set variables to control iterative procedure
  Target<-length(Y) - ncol(X) # Target value of iterations. RHS of Eq 21.
  GInt <- 30 # Number of intervals to consider (doubled below)
  Gstep <- 2*var(Y)/(GInt-1) # Length of step
  ## Coarse intervals.
  Omega <- diag(nrow(X)) # OLS weighting matrix (identity), just to initialize size
  Iterations <- matrix(0,ncol=2,nrow=GInt) # Empty dataframe to store results
  Iterations <- data.frame(Iterations); names(Iterations) <- c('GSQ','Deviation') # Formatting empty dataframe.  Will contain estimated model-error variance (GSQ) and deviation from target.
  for (w in 1:GInt) {
    GSQ <- (w-1)*Gstep; Iterations$GSQ[w]<-GSQ; # guess at model error variance
    for (i in 1:length(Y)) {
      for (j in 1:length(Y)) {
        if (i==j) { # Diagonal elements of weighting matrix
          if (SkewAdj) { # if skew-adjusted GLS, use Eq 22 (correct manual to match Giffis and Stedinger (2007) Eq 5; WREG v1.05 code is correct)
            Omega[i,j] <-  GSQ + 
              Yhat.SigReg[i]^2*(1+LP3$K[i]*LP3$G[i]+
                  0.5*LP3$K[i]^2*(1+0.75*LP3$G[i]^2)+
                  Wg[i]*LP3$K[i]*dKdG[i]*(3*LP3$G[i]+0.75*LP3$G[i]^3)+
                  Wg[i]^2*dKdG[i]^2*(6+9*LP3$G[i]^2+1.875*LP3$G[i]^4))/M[i,j] + 
              (1-Wg[i])^2*Yhat.SigReg[i]^2*MSEGR*dKdG[i]^2
          } else { # if normal GLS, use Eq 19.
            Omega[i,j] <-  GSQ + Yhat.SigReg[i]^2*
              (1+LP3$K[i]*LP3$G[i]+0.5*LP3$K[i]^2*
                  (1+0.75*LP3$G[i]^2))/M[i,j]
          }
        } else { # Off-diagonal elements of weighting matrix
          if (SkewAdj) { # if skew-adjusted GLS, use Eq 22 (correct manual to match Giffis and Stedinger (2007) Eq 5; WREG v1.05 code is correct)
            Omega[i,j] <- Rhos[i,j]*Yhat.SigReg[i]*Yhat.SigReg[j]*M[i,j]*
              (1+0.5*LP3$K[i]*LP3$G[i]+0.5*LP3$K[j]*LP3$G[j]+
                  0.5*LP3$K[i]*LP3$K[j]*(Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                  0.5*Wg[j]*LP3$K[i]*LP3$G[j]*dKdG[j]*(3*Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                  0.5*Wg[i]*LP3$K[j]*LP3$G[i]*dKdG[i]*(3*Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                  Wg[i]*Wg[j]*Yhat.SigReg[i]*Yhat.SigReg[j]*dKdG[i]*dKdG[j]*Covgg[i,j])/M[i,i]/M[j,j]
          } else { # if normal GLS, use Eq 19.
            Omega[i,j] <- Rhos[i,j]*Yhat.SigReg[i]*Yhat.SigReg[j]*M[i,j]*
              (1+0.5*LP3$K[i]*LP3$G[i]+0.5*LP3$K[j]*LP3$G[j]+
                  0.5*LP3$K[i]*LP3$K[j]*(Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j]))/M[i,i]/M[j,j]
          }
        }
      }
    }
    B_hat <- solve(t(X)%*%solve(Omega)%*%X)%*%t(X)%*%solve(Omega)%*%Y # Use weighting matrix to estimate regression coefficients
    Iterations$Deviation[w] <- t(Y-X%*%B_hat)%*%solve(Omega)%*%(Y-X%*%B_hat)-Target # Difference between result and target value. Eq 21.
  }
  ## Finer intervals.  Expands iterval with sign change to get closer to zero.
  Signs <- sign(Iterations$Deviation); LastPos <- which(diff(Signs)!=0) # Finds where sign changes from positive to negative
  if (length(LastPos)==0) { # If no change in sign of deviation
    BestPos <- which(abs(Iterations$Deviation)==min(abs(Iterations$Deviation))) # Use the minimum deviation
    GSQ <- Iterations$GSQ[BestPos] # Best estiamte of model-error variance
  } else { # There is a sign change, so expand the interval with sign change.
    Iterations2 <- matrix(0,ncol=2,nrow=GInt) # Empty matrix to store finer intervals.
    Iterations2 <- data.frame(Iterations); names(Iterations) <- c('GSQ','Deviation') # Formatting empty matrix
    Gstep <- (Iterations$GSQ[LastPos+1]-Iterations$GSQ[LastPos])/(GInt-1) # Step for finer intervals
    for (w in 1:GInt) {
      GSQ <- Iterations$GSQ[LastPos]+(w-1)*Gstep; Iterations2$GSQ[w]<-GSQ; # estimate of model-error variance
      for (i in 1:length(Y)) {
        for (j in 1:length(Y)) {
          if (i==j) { # Diagonal elements of weighting matrix
            if (SkewAdj) { # if skew-adjusted GLS, use Eq 22 (correct manual to match Giffis and Stedinger (2007) Eq 5; WREG v1.05 code is correct)
              Omega[i,j] <-  GSQ + 
                Yhat.SigReg[i]^2*(1+LP3$K[i]*LP3$G[i]+
                    0.5*LP3$K[i]^2*(1+0.75*LP3$G[i]^2)+
                    Wg[i]*LP3$K[i]*dKdG[i]*(3*LP3$G[i]+0.75*LP3$G[i]^3)+
                    Wg[i]^2*dKdG[i]^2*(6+9*LP3$G[i]^2+1.875*LP3$G[i]^4))/M[i,j] + 
                (1-Wg[i])^2*Yhat.SigReg[i]^2*MSEGR*dKdG[i]^2
            } else { # if normal GLS, use Eq 19.
              Omega[i,j] <-  GSQ + Yhat.SigReg[i]^2*
                (1+LP3$K[i]*LP3$G[i]+0.5*LP3$K[i]^2*
                    (1+0.75*LP3$G[i]^2))/M[i,j]
            }
          } else { # Off-diagonal elements of weighting matrix
            if (SkewAdj) { # if skew-adjusted GLS, use Eq 22 (correct manual to match Giffis and Stedinger (2007) Eq 5; WREG v1.05 code is correct)
              Omega[i,j] <- Rhos[i,j]*Yhat.SigReg[i]*Yhat.SigReg[j]*M[i,j]*
                (1+0.5*LP3$K[i]*LP3$G[i]+0.5*LP3$K[j]*LP3$G[j]+
                    0.5*LP3$K[i]*LP3$K[j]*(Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                    0.5*Wg[j]*LP3$K[i]*LP3$G[j]*dKdG[j]*(3*Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                    0.5*Wg[i]*LP3$K[j]*LP3$G[i]*dKdG[i]*(3*Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                    Wg[i]*Wg[j]*Yhat.SigReg[i]*Yhat.SigReg[j]*dKdG[i]*dKdG[j]*Covgg[i,j])/M[i,i]/M[j,j]
            } else { # if normal GLS, use Eq 19.
              Omega[i,j] <- Rhos[i,j]*Yhat.SigReg[i]*Yhat.SigReg[j]*M[i,j]*
                (1+0.5*LP3$K[i]*LP3$G[i]+0.5*LP3$K[j]*LP3$G[j]+
                    0.5*LP3$K[i]*LP3$K[j]*(Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j]))/M[i,i]/M[j,j]
            }
          }
        }
      }
      B_hat <- solve(t(X)%*%solve(Omega)%*%X)%*%t(X)%*%solve(Omega)%*%Y # Use weighting matrix to estimate regression coefficients
      Iterations2$Deviation[w] <- t(Y-X%*%B_hat)%*%solve(Omega)%*%(Y-X%*%B_hat)-Target # Difference between result and target value. Eq 21.
    }
    BestPos <- which(Iterations2$Deviation==min(Iterations2$Deviation[Iterations2$Deviation>0])) # Find minimum positive deviation from Eq 21.
    GSQ <- Iterations2$GSQ[BestPos] # best estimate of model-error variance
  }
  ## Calculate Final Omega (weighting matrix)
  for (i in 1:length(Y)) {
    for (j in 1:length(Y)) {
      if (i==j) { # Diagonal elements of weighting matrix
        if (SkewAdj) { # if skew-adjusted GLS, use Eq 22 (correct manual to match Giffis and Stedinger (2007) Eq 5; WREG v1.05 code is correct)
          Omega[i,j] <-  GSQ + 
            Yhat.SigReg[i]^2*(1+LP3$K[i]*LP3$G[i]+
                0.5*LP3$K[i]^2*(1+0.75*LP3$G[i]^2)+
                Wg[i]*LP3$K[i]*dKdG[i]*(3*LP3$G[i]+0.75*LP3$G[i]^3)+
                Wg[i]^2*dKdG[i]^2*(6+9*LP3$G[i]^2+1.875*LP3$G[i]^4))/M[i,j] + 
            (1-Wg[i])^2*Yhat.SigReg[i]^2*MSEGR*dKdG[i]^2
        } else { # if normal GLS, use Eq 19.
          Omega[i,j] <-  GSQ + Yhat.SigReg[i]^2*
            (1+LP3$K[i]*LP3$G[i]+0.5*LP3$K[i]^2*
                (1+0.75*LP3$G[i]^2))/M[i,j]
        }
      } else { # Off-diagonal elements of weighting matrix
        if (SkewAdj) { # if skew-adjusted GLS, use Eq 22 (correct manual to match Giffis and Stedinger (2007) Eq 5; WREG v1.05 code is correct)
          Omega[i,j] <- Rhos[i,j]*Yhat.SigReg[i]*Yhat.SigReg[j]*M[i,j]*
            (1+0.5*LP3$K[i]*LP3$G[i]+0.5*LP3$K[j]*LP3$G[j]+
                0.5*LP3$K[i]*LP3$K[j]*(Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                0.5*Wg[j]*LP3$K[i]*LP3$G[j]*dKdG[j]*(3*Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                0.5*Wg[i]*LP3$K[j]*LP3$G[i]*dKdG[i]*(3*Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j])+
                Wg[i]*Wg[j]*Yhat.SigReg[i]*Yhat.SigReg[j]*dKdG[i]*dKdG[j]*Covgg[i,j])/M[i,i]/M[j,j]
        } else { # if normal GLS, use Eq 19.
          Omega[i,j] <- Rhos[i,j]*Yhat.SigReg[i]*Yhat.SigReg[j]*M[i,j]*
            (1+0.5*LP3$K[i]*LP3$G[i]+0.5*LP3$K[j]*LP3$G[j]+
                0.5*LP3$K[i]*LP3$K[j]*(Rhos[i,j]+0.75*LP3$G[i]*LP3$G[j]))/M[i,i]/M[j,j]
        }
      }
    }
  }
  ## Control output
  GLS.Weights <- list(GSQ=GSQ,Omega=Omega) # Output contains model-error variance and weighting matrix.
  return(GLS.Weights)
}


#' Exploring correlation function
#'
#'@description
#'\code{testXCorr} explores the fit of the proposed approximation of 
#'cross-correlation.
#'
#' @param siteData A list of the records at each site.  Each element of the list,
#' representing a site, is a data frame with 2 variables names \code{year} 
#' and \code{obs}.
#' @param BasinChars A dataframe containing three variables: \code{StationID} is the 
#' numerical identifier (without a leading zero) of each site, \code{Lat} is the latitude
#' of the site, in decimal degrees, and \code{Long} is the longitude of the site, in decimal
#' degrees.
#' @param alpha A number, \code{alpha} is a parameter used in the estimated 
#' cross-correlation between site records.  See equation 20 in the WREG 
#' v. 1.05 manual.
#' @param theta A number, \code{theta} is a parameter used in the estimated 
#' cross-correlation between site records.  See equation 20 in the WREG 
#' v. 1.05 manual.
#' @param DistMeth A value of \code{1} indicates that the "Nautical Mile" 
#' approximation should be used to calculate inter-site distances.  A value 
#' of \code{2} designates the Haversine approximation.  See 
#' \code{\link{Dist.WREG}}.
#' @param concurrentMin A number specifying the minimum number of years of 
#' concurrent record required to estimate cross-correlation.
#' @param plot A logical spcifying if the plot should be created.
#' 
#' @return If \code{plot=FALSE}, the Nash-Sutcliffe model efficieny is returned.
#' 
#' @export
testXCorr <- function(siteData,BasinChars,alpha,theta,DistMeth=c(1,2),
  concurrentMin=25,plot=TRUE) {
  # William Farmer, October 22, 2015
  #
  
  plotData <- matrix(NA,ncol=2,nrow=length(siteData)^2)
  RecordLengths <- matrix(NA,ncol=length(siteData),nrow=length(siteData))
  iter <- 0
  maxDist <- maxcor <- -Inf
  mincor <- Inf
  for (i in 1:length(siteData)) {
    iData <- siteData[[i]]
    for (j in i:length(siteData)) {
      if (i!=j) {
        ijDist <- Dist.WREG(BasinChars$Lat[i],BasinChars$Long[i],
          BasinChars$Lat[j],BasinChars$Long[j],method=DistMeth)
        maxDist <- max(maxDist,ijDist)
        jData <- siteData[[j]]
        ### Concurrent period of record
        interY <- intersect(iData$year,jData$year)
        RecordLengths[i,j] <- length(interY)
        if (RecordLengths[i,j]>2) {
          indx <- which(is.element(iData$year,interY))
          jndx <- which(is.element(jData$year,interY))
          ijcor <- NA
          if (sd(iData$obs[indx])>0&sd(jData$obs[jndx])>0) {
            ijcor <- cor(iData$obs[indx],jData$obs[jndx],
              method='pearson')
          }
          mincor <- min(mincor,ijcor,na.rm=T)
          maxcor <- max(maxcor,ijcor,na.rm=T)
          if (RecordLengths[i,j]>=concurrentMin) {
            # if above minimum, add to plot data
            iter <- iter + 1
            plotData[iter,1] <- ijDist
            plotData[iter,2] <- ijcor
          }
        }
      }
    }
  }
  
  ndx <- which(!is.na(rowSums(plotData)))
  plotData <- plotData[ndx,]
  
  estRhos <- theta^(plotData[,1]/(alpha*plotData[,1]+1))
  
  nse <- 1 - sum((estRhos-plotData[,2])^2)/
    sum((plotData[,2]-mean(plotData[,2]))^2)
  
  if (plot) {
    ylim1 <- ifelse(round(mincor,digits=1)>mincor,round(mincor,digits=1)-0.1,
      round(mincor,digits=1))
    ylim2 <- ifelse(round(maxcor,digits=1)<maxcor,round(maxcor,digits=1)+0.1,
      round(maxcor,digits=1))
    ylim <- c(ylim1,ylim2)
    ny <- round((ylim2-ylim1)/0.1)
    nx <- 10
    splits <- c(1,5,10,50,100,500,1000,5000,10000)
    splits1 <- splits-(maxDist/nx)
    splits1[splits1<0] <- NA
    ndx <- which(splits1==min(splits1,na.rm=T))
    xlim2 <- splits[ndx]*nx
    plotDists <- seq(0,xlim2,length.out=1000)
    plotRhos <- theta^(plotDists/(alpha*plotDists+1))
    xlim <- c(0,xlim2)
    plot(plotData[,1],plotData[,2],type='p',
      main=paste0('Correlation Smoothing Function\n',
        '(alpha=',alpha,', theta=',theta,', NSE=',round(nse,digits=4),')'),
      xlab='Geographic Distance (km)',
      ylab='Sample rho',
      xaxs='i',yaxs='i',
      ylim=ylim,xlim=xlim)
    lines(plotDists,plotRhos,lty=5,col='red')
    grid(nx=nx,ny=ny)
    legend('topright',
      legend=c(paste0('Observation (Y=',concurrentMin,')'),'Model'),lty=c(NA,5),
      col=c('black','red'),pch=c(0,NA))
  } else {
    return(nse)
  }
  
  
  
}