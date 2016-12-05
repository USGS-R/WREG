#'Leverage Statistics (WREG)
#'
#'@description The \code{Leverage} function calculates the leverage statistics
#'for each observation in the regression.
#'
#'@param X contains independent variables.  A leading constant of one is
#'  included if a constant is included in the regression.  Each row represents a
#'  unique observation.
#'@param Omega is the weighting matrix used for regression fitting
#'@param Ch allows the user to specify a custom criteria for testing.  If not
#'  specified, the default is 4 for region-of-influence regression and 2
#'  elsewise.
#'@param x0 is used only in the case of region-of-influence regression.  It
#'  contains the independent variables of a specific site against which to
#'  calculate the leverage.
#'@param ROI is a logical vector specifying if the regression is or is not 
#'  region-of-influence regression
#'  
#'@details Leverage is a measure of how far way the independent variables of an
#'observation are from the other observations in the regression set.  A leverage
#'is considered significant if the absolute value of the leverage is greater
#'than the critical value.  These calculations are based on equations 40, 41 and
#'42 of the WREG v. 1.0 manual.
#'
#'@return The function returns a list as output.  This list contains: 
#'  \item{Leverage}{A vector indicating the leverage of each observation. This
#'  is a vector whose length equals the number of rows in X.} \item{Limit}{A
#'  number indicating the critical value of leverage for this data set.} 
#'  \item{Significant}{A logical vector the same size as \code{Leverage}. It
#'  indicates if the leverage is significant for each observation.}
#'  
#'@examples
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
#' # calculate leverage of each point
#' leverageResult <- Leverage(X = X, 
#'   Omega = result$Weighting)
#'
#'@export
Leverage <- function(X,Omega,Ch=NA,x0=NA,ROI=FALSE) {
  # William Farmer, USGS, January 02, 2015
  # 01/27/2015, WHF: Added ABS on limits
  # 11/10/16 Greg Petrochenkov: Changed validation scheme
  
  # Some upfront error handling
  if(!wregValidation(missing(X), "eq", FALSE,
                     "Independent variables (X) must be provided.", warnFlag = TRUE)){
  
    if(!wregValidation(X, "numeric", message="X must be numeric", warnFlag = TRUE)){
      
      if(!wregValidation((length(unique(apply(X,FUN=class,MARGIN=2)))!=1)|
                         (unique(apply(X,FUN=class,MARGIN=2))!="numeric"), "eq", FALSE,
                         "Independent variables (X) must be provided as class numeric.", warnFlag = TRUE)){
        
        if(!wregValidation(sum(is.na(as.matrix(X))), "eq", 0,
                           paste0("Some independent variables (X) contain missing ",
                                  "values.  These must be removed."), warnFlag = TRUE)){
          
          wregValidation(sum(is.infinite(as.matrix(X))), "eq", 0,
                         paste0("Some independent variables (X) contain infinite ",
                                "values.  These must be removed."), warnFlag = TRUE)
        }
      }
    }
  }
  if(!wregValidation(missing(Omega), "eq", FALSE,
                     "A weighting matrix (Omega) must be provided.", warnFlag = TRUE)){
    
    if(!wregValidation(!is.matrix(Omega), "eq", FALSE,
                       "The weighting matrix (Omega) must be provided as a matrix.", warnFlag = TRUE)){
      
      if(!wregValidation(Omega, "numeric", message =
                         "The weighting matrix (Omega) must be provided as class numeric.", warnFlag = TRUE)){
        
        wregValidation(det(Omega), "notEq", 0,
                       paste0("Some independent variables (X) contain infinite ",
                              "values.  These must be removed."), warnFlag = TRUE)
      }
    }
  }
  
  
  if (warn("check")) {
    stop("Invalid inputs were provided.  See warnings().", warn("get"))
  }
  
  if (ROI==F) { # Not region-of-influence regression
    h <- diag(as.matrix(X)%*%solve(t(X)%*%solve(Omega)%*%as.matrix(X))%*%t(X)%*%solve(Omega)) # Leverage, Eq 40
    if (is.na(Ch)) {Ch=2} # Default for non-ROI 
  } else if (ROI==T) { # Region-of-influence regression
    h <- t(as.matrix(x0)%*%solve(t(X)%*%solve(Omega)%*%as.matrix(X))%*%t(X)%*%solve(Omega)) # Leverage, Eq 41
    if (is.na(Ch)) {Ch=4} # Default for ROI 
  }
  h_limit <- Ch*mean(h) # Critical value of leverage, Eq 42
  sig <- abs(h)>abs(h_limit) # Significance test. Logical: Leverage is signifianct.
  Output <- list(Leverage=h,Limit=h_limit,Significant=sig)
  return(Output)
}