#'Weighted-Multiple-Linear Regression Program (WREG)
#'
#'@description The \code{WREG.UW} function executes the multiple linear 
#'  regression analysis using a user-provided weighting matrix.
#'  
#'@param Y The dependent variable of interest, with any transformations already 
#'  applied.
#'@param X The independent variables in the regression, with any transformations
#'  already applied.  Each row represents a site and each column represents a 
#'  particular independe variable.  (If a leading constant is used, it should be
#'  included here as a leading column of ones.)  The rows must be in the same 
#'  order as the dependent variables in \code{Y}.
#'@param customWeight This allows the user to enter a custom weighting matrix. 
#'  It is included also to provide legacy code for WREG v. 1.05. 
#'  \code{customWeight} can either be a square matrix of weights with a length 
#'  equal to \code{length(Y)} or a list containing three elements.  The elements
#'  of the list include (1) \code{Omega} as the square weighting matrix, (2) 
#'  \code{var.modelerror.k} as the estimated variance of the model errors from 
#'  the k-variable model (\code{k=ncol(X)}), and (3) \code{var.modelerror.0} as 
#'  the variance of the model errors from a consant-only regression.
#'@param transY A required character string indicating if the the 
#'  dependentvariable was transformed by the common logarithm ('log10'), 
#'  transformed by the natural logarithm ('ln') or untransformed ('none').
#'@param x0 A vector containing the independent variables (as above) for a 
#'  particular target site.  This variable is only used for ROI analysis.
#'  
#'  
#'@details This function allows users to develop weights outside of the WREG 
#'  program and observe the resultant regressions.  Note that the weighting
#'  matrix must be invertible.
#'  
#'@return All outputs are returned as part of a list.  The elements of the list 
#'  depend on the type of regression performed.  The elements of the list may 
#'  include: \item{Coefs}{A data frame composed of four variables: (1) 
#'  \code{Coefficient} contains the regression coefficeints estimated for the 
#'  model, (2) \code{Standard Error} contains the standard errors of 
#'  each regression coefficient, (3) \code{tStatistic} contains the Student's 
#'  T-statistic of each regression coefficient and (4) \code{pValue} contains 
#'  the significance probability of each regression coefficient.} 
#'  \item{ResLevInf}{A data frame composed of three variables for each site in 
#'  the regression.  \code{Residual} contains the model residuals. 
#'  \code{Leverage} contains the leverage of each site.  \code{Influence} 
#'  contains the influence of each site.} \item{LevLim}{The critical value of 
#'  leverage.  See \code{\link{Leverage}}} \item{InflLim}{The critical value of 
#'  influence.  See \code{\link{Influence}}} \item{LevInf.Sig}{A logical matrix 
#'  indicating if the leverage (column 1) is significant and the influence 
#'  (column 2) is significant for each site in the regression.} 
#'  \item{PerformanceMetrics}{A list of not more than ten elements.  All 
#'  regression types return the mean squared error of residuals (\code{MSE}), 
#'  the coefficient of determination (\code{R2}), the adjusted coefficient of 
#'  determination (\code{R2_adj}) and the root mean squared error (\code{RMSE}, 
#'  in percent).  If \code{customWeight} contains model error variances, then
#'  the pseudo coefficient of regression (\code{R2_pseudo}), the average
#'  variance of prediction (\code{AVP}), the standard error of prediction
#'  (\code{Sp}, in percent), a vector of the individual variances of prediction
#'  for each site (\code{VP.PredVar}), the model-error variance
#'  (\code{ModErrVar}) and the standardized model error variance
#'  (\code{StanModErr}, in percent) are also returned.  Details on the
#'  appropriateness and applicability of performance metrics can be found in the
#'  WREG manual.} \item{X}{The input predictors.} \item{Y}{The input 
#'  observations.} \item{fitted.values}{A vector of model estimates from the 
#'  regression model.} \item{residuals}{A vector of model residuals.} 
#'  \item{Weighting}{The weighting matrix used to develop regression estimates.}
#'  \item{Input}{A list of input parameters for error searching.  Currently 
#'  empty.}
#'@import stats
#'  
#' @examples
#' # Import some example data
#' peakFQdir <- paste0(
#'   file.path(system.file("exampleDirectory", package = "WREG"),
#'     "pfqImport"))
#' gisFilePath <- file.path(peakFQdir, "pfqSiteInfo.txt")
#' importedData <- importPeakFQ(pfqPath = peakFQdir, gisFile = gisFilePath)
#' 
#' # Organizing input data
#' Y <- importedData$Y$AEP_0.5
#' X <- importedData$X[c("Sand", "OutletElev", "Slope")]
#' transY <- "none"
#' 
#' # Make simple weighting using inverse record lengths
#' inverseRecLen <- diag(1 / diag(importedData$recLen))
#' 
#' # Run user-weights regression
#' result <- WREG.UW(Y, X, customWeight = inverseRecLen, transY)
#' 
#'@export
WREG.UW <- function(Y,X,customWeight,transY,x0=NA) {
  # William Farmer, USGS, January 05, 2015
  # Greg PEtrochenkov, USGS, November 14, 2016 :  Changed validation scheme

  warn("clear")
  # Some upfront error handling
  wregValidation((!missing(X)&!missing(Y))&&(length(Y)!=nrow(X)), "eq", FALSE,
                 paste0("The length of Y must be the same as ",
                        "the number of rows in X."), warnFlag = TRUE)
  
  if (!wregValidation(missing(Y), "eq", FALSE,
                      "Dependent variable (Y) must be provided", warnFlag = TRUE)) {
    
    if (!wregValidation(Y, "numeric", message = 
                        "Dependent variable (Y) must be provided as class numeric",
                        warnFlag = TRUE)) {
      
      wregValidation(sum(is.na(Y)), "eq", 0 ,
                     paste0("The depedent variable (Y) contains missing ",
                            "values.  These must be removed."),
                     warnFlag = TRUE)
      
      wregValidation(sum(is.infinite(Y)), "eq", 0 ,
                     paste0("The depedent variable (Y) contains infinite ",
                            "values.  These must be removed."),
                     warnFlag = TRUE)
    }
  }
  
  if (!wregValidation(missing(X), "eq", FALSE,
                      "Independent variables (X) must be provided.", warnFlag = TRUE)) {
    
    if (!wregValidation((length(unique(apply(X,FUN=class,MARGIN=2)))!=1)|
                        (unique(apply(X,FUN=class,MARGIN=2))!="numeric"), "eq", FALSE,
                        "Independent variables (X) must be provided as class numeric.", warnFlag = TRUE)){
      
      wregValidation(sum(is.na(as.matrix(X))), "eq", 0,
                     paste0("Some independent variables (X) contain missing ",
                            "values.  These must be removed."), warnFlag = TRUE)
      
      wregValidation(sum(is.infinite(as.matrix(X))), "eq", 0,
                     paste0("Some independent variables (X) contain infinite ",
                            "values.  These must be removed."), warnFlag = TRUE)
    }
  }
  
  if(!wregValidation(missing(transY)|!is.character(transY),"eq", FALSE,
                     "transY must be included as a character string.", warnFlag=TRUE)){
    
    wregValidation(!is.element(transY,c("none","log10","ln")), "eq", FALSE,
                   "transY must be either 'none', 'log10' or 'ln'.", warnFlag = TRUE)
  }
  
  wregValidation(missing(customWeight)|(!is.matrix(customWeight)&!is.list(customWeight)), "eq", FALSE,
                 "Custom weighting matrix must be provided as a list or matrix.", warnFlag = TRUE)
  
  ## Determine if ROI is being applied
  if (is.na(sum(x0))) { # ROI regression is not used.
    ROI <- F
  } else { # ROI regression is used.
    ROI <- T
  }
  ## Just initial values for control.
  var.modelerror.k <- NA
  customModelError <- FALSE
  
  ## Weighting matrix
  # Allows user to specify particular weighting scheme.  Useful for legacy code.
  if (is.list(customWeight)) { # Omega from legacy code; Also contains information on model-error variances.
    Omega <- customWeight$Omega # Custom weighting
    var.modelerror.k <- customWeight$var.modelerror.k # Custom k-variable model-error variance
    var.modelerror.0 <- customWeight$var.modelerror.0 # Custom constant-model model-error variance
    customModelError=TRUE # Logical to note that the user has provided information on the model-error variance
  } else { # User-defined weighting, with no informaiton on model-error variance.
    Omega <- customWeight # Custom weighting
    var.modelerror.k <- NA # NULL custom k-variable model-error variance to control for errors.
    var.modelerror.0 <- NA # NULL custom constant-model model-error variance to control for errors.
  }
  
  wregValidation(det(Omega), "notEq", 0,
                 paste("The weighting matrix is singular and, therefore,",
                       "cannot be inverted.  Reconsider the weighting matrix."), warnFlag = TRUE)
  
  if (warn("check")) {
    stop('Invalid inputs were provided. See warnings().', warn("get"))
  }
  
  #Convert X and Y from dataframes to matrices to work with matrix operations below
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
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
  RMSE <- NA
  if (transY=='log10') {
    RMSE <-100*sqrt(exp(log(10)*log(10)*MSE)-1) # Root-mean-squared error, in percent. Eq 34
  } else if (transY=='ln') {
    RMSE <-100*sqrt(exp(MSE)-1) # Root-mean-squared error, in percent. transformed for natural logs.
  }
  PerfMet <- list(MSE=MSE,R2=R2,R2_adj=R2_adj,RMSE=RMSE) # Performance metrics for output (basic, for OLS)
  if (customModelError==T) { # non-OLS requires additional performance metrics
    R2_pseudo <- 1 - var.modelerror.k/var.modelerror.0 # Pseudo coefficient of determination. Eq 39
    AVP <- var.modelerror.k + mean(diag(X%*%solve(t(X)%*%solve(Omega)%*%X)%*%t(X))) # Average varaince of prediction. Eq 32
    VP <- vector(length=length(Y)) # Empty vector for individual variances of prediction
    for (i in 1:length(VP)) { # Individual variances of prediction
      VP[i] <- var.modelerror.k + X[i,]%*%solve(t(X)%*%solve(Omega)%*%X)%*%X[i,] # Individual variance of prediction.  Based on Eq 32.
    }
    VP <- data.frame(VP); names(VP) <- 'PredVar' # Formating the VP vector for output
    Sp <- Se <- NA
    if (transY=='log10') {
      Sp <- 100*sqrt(exp(log(10)*log(10)*AVP)-1) # Standard error of predictions. Eq 33.
      Se <-100*sqrt(exp(log(10)*log(10)*var.modelerror.k)-1) # Standard model error. Not noted in the manual, but included as output in WREG 1.05.  Based on Eq 33.
    } else if (transY=='ln') {
      # corrected for natural logs
      Sp <- 100*sqrt(exp(AVP)-1)
      Se <-100*sqrt(exp(var.modelerror.k)-1)
    }
    PerfMet <- c(PerfMet,R2_pseudo=R2_pseudo,AVP=AVP,Sp=Sp,VP=VP,ModErrVar=var.modelerror.k,StanModErr=Se) # Performance metrics for output
  }
  
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
  Output <- list(Coefs=Coefs,ResLevInf=ResLevInf,LevLim=Lev$Limit,
    InflLim=Infl$Limit,LevInf.Sig=LevInf.Sig,
    PerformanceMetrics=PerfMet,X=X,Y=Y,fitted.values=Y_hat,residuals=e,
    Weighting=Omega,Inputs=list(transY=transY))
  if (ROI) { # Appended at-site estimates for ROI calculations
    Y_est <- as.matrix(x0)%*%B_hat # ROI site estimate
    Output <- c(Output,Y.ROI=Y_est,x0.ROI=x0)
  }
  
  class(Output) <- 'WREG.UW'
  
  return(Output)
}