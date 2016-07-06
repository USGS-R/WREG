#'Weighted-Multiple-Linear Regression Program (WREG)
#'
#'@description The \code{WREG.OLS} function executes the multiple linear 
#'  regression analysis using ordinary least-squares regression.
#'  
#'@param Y The dependent variable of interest, with any transformations already 
#'  applied.
#'@param X The independent variables in the regression, with any transformations
#'  already applied.  Each row represents a site and each column represents a 
#'  particular independe variable.  (If a leading constant is used, it should be
#'  included here as a leading column of ones.)  The rows must be in the same 
#'  order as the dependent variables in \code{Y}.
#'@param transY A required character string indicating if the the 
#'  dependentvariable was transformed by the common logarithm ('log10'), 
#'  transformed by the natural logarithm ('ln') or untransformed ('none').
#'@param x0 A vector containing the independent variables (as above) for a 
#'  particular target site.  This variable is only used for ROI analysis.
#'  
#'@details This function follows the basic implementation of ordinary 
#'  least-squares regression.
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
#'  in percent).  Details on the appropriateness and applicability of
#'  performance metrics can be found in the WREG manual.} \item{X}{The input
#'  predictors.} \item{Y}{The input observations.} \item{fitted.values}{A vector
#'  of model estimates from the regression model.} \item{residuals}{A vector of
#'  model residuals.} \item{Weighting}{The weighting matrix used to develop
#'  regression estimates.} \item{Input}{A list of input parameters for error
#'  searching.  Currently empty.}
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
#' # Run a simple regression
#' Y <- importedData$Y$AEP_0.5
#' X <- importedData$X[c("Sand", "OutletElev", "Slope")]
#' transY <- "none"
#' result <- WREG.OLS(Y, X, transY)
#'
#'@export

WREG.OLS <- function(Y,X,transY,x0=NA) {
  # William Farmer, USGS, January 05, 2015
  
  # Some upfront error handling
  err <- FALSE
  if ((!missing(X)&!missing(Y))&&
      (length(Y)!=nrow(X))) {
    warning(paste0("The length of Y must be the same as ",
      "the number of rows in X."))
    err <- TRUE
  }
  if (missing(Y)) {
    warning("Dependent variable (Y) must be provided.")
    err <- TRUE
  } else {
    if (!is.numeric(Y)) {
      warning("Dependent variable (Y) must be provided as class numeric.")
      err <- TRUE
    } else {
      if (sum(is.na(Y))>0) {
        warning(paste0("The depedent variable (Y) contains missing ",
          "values.  These must be removed."))
        err <- TRUE
      }
      if (sum(is.infinite(Y))>0) {
        warning(paste0("The depedent variable (Y) contains infinite ",
          "values.  These must be removed."))
        err <- TRUE
      }
    }
  }
  if (missing(X)) {
    warning("Independent variables (X) must be provided.")
    err <- TRUE
  } else {
    if ((length(unique(apply(X,FUN=class,MARGIN=2)))!=1)|
        (unique(apply(X,FUN=class,MARGIN=2))!="numeric")) {
      warning("Independent variables (X) must be provided as class numeric.")
      err <- TRUE
    } else {
      if (sum(is.na(as.matrix(X)))>0) {
        warning(paste0("Some independent variables (X) contain missing ",
          "values.  These must be removed."))
        err <- TRUE
      }
      if (sum(is.infinite(as.matrix(X)))>0) {
        warning(paste0("Some independent variables (X) contain infinite ",
          "values.  These must be removed."))
        err <- TRUE
      }
    }
  }
  if(missing(transY)|!is.character(transY)) {
    warning("transY must be included as a character string.")
    err <- TRUE
  } else if (!is.element(transY,c("none","log10","ln"))) {
    warning("transY must be either 'none', 'log10' or 'ln'.")
    err <- TRUE
  }

  ## Determine if ROI is being applied
  if (is.na(sum(x0))) { # ROI regression is not used.
    ROI <- F
  } else { # ROI regression is used.
    ROI <- T
  }
  if (ROI&&(length(x0)!=ncol(X))) {
    warning(paste0("The length of x0 must be the same as ",
      "the number of columns in X"))
    err <- TRUE
  }
  if (ROI&&(!is.numeric(x0))) {
    warning(paste0("The input x0 must be of the numeric class"))
    err <- TRUE
  }
  if (err) {
    stop('Invalid inputs were provided. See warnings().')
  }
  ## Just initial values for control.
  var.modelerror.k <- NA

  ## Weighting matrix
  # Ordinary Least Squares
  Omega <- diag(nrow(X)) # weighting matrix
  
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

  ## Leverage and influence statistics
  Lev <- Leverage(X=X,Omega=Omega,x0=x0,ROI=ROI) # Leverage subroutine
  Infl <- Influence(e=e,X=X,Omega=Omega,Beta=B_hat,ROI=ROI,Lev=Lev$Leverage) # Influence subroutine
  
  ## Significance of regression parameters
  B_var <- diag(solve(t(X)%*%solve(Omega)%*%X)) # Covariances of regression coefficients. Eq 46
  # Eq 46 in Manual for v1.05 is incoorect.  As reflected in code for v1.05 and independent verification, formula is altered for OLS.
    B_var <- MSE*B_var # Altered Eq 46 for OLS
  
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
    InflLim=Infl$Limit,LevInf.Sig=cbind(Lev$Significant,Infl$Significant),
    PerformanceMetrics=PerfMet,X=X,Y=Y,fitted.values=Y_hat,residuals=e,
    Weighting=Omega,Inputs=list(transY=transY))
  if (ROI) { # Appended at-site estimates for ROI calculations
    Y_est <- as.matrix(x0)%*%B_hat # ROI site estimate
    Output <- c(Output,Y.ROI=Y_est,x0.ROI=x0)
  }
  
    class(Output) <- 'WREG.OLS'
    
    return(Output)
}