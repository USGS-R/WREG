#' @title print.WREG
#' @description Print methods for WREG output lists
#' @name print
#' @aliases print.OLS
#' @aliases print.WLS
#' @aliases print.GLS
#' @aliases print.GLSs
#'
#' @title print.WREG
#'
#' @param x An output list from one of the WREG... functions
#'
#' @note \code{print} is a generic name for the functions documented.
#' \cr
#' If called, \code{print} displays a summary of the output from WREG... functions
#'
#' @rdname print
#' @export
#' 
#' @rdname print.OLS
#' @return \code{print.WREG.OLS} Prints a summary of output list from WREG.OLS
#' @examples 
#' ## print.WREG.OLS
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
#' print(result)
#' @export
print.WREG.OLS <- function(x, ...) {
  object <- x
  cat(paste0("Regression Model for ",names(object$Y),'\n'))
  cat('Coefficients fit by ordinary least-squares.\n')
  cat("\nCall:\n", object$Inputs$call,'\n\n')
  cat(paste0('\nPerformance Metrics\n(Note: Units are based',
             ' on the transformation of the dependent variable.)\n'))
  cat('Mean Squared Error:\t',
      round(object$PerformanceMetrics$MSE,digits=4),'\n')
  cat(paste0('Root Mean Squared Error:\t',
             round(object$PerformanceMetrics$RMSE,digits=2),'%\n'))
  cat(paste0('Coefficient of Determination (R2):\t',
             round(100*object$PerformanceMetrics$R2,digits=2),'%\n'))
  cat(paste0('Adjusted Coefficient of Determination (Adj-R2):\t',
             round(100*object$PerformanceMetrics$R2_adj,digits=2),'%\n\n'))
  cat('Coefficients of Model')
  object$Coefs[,1:3] <- round(object$Coefs[,1:3],3)
  object$Coefs[,4] <- round(object$Coefs[,4],4)
  object$Coefs[,5] <- ifelse(object$Coefs[,4]<=0.05,'*','')
  names(object$Coefs) <- c('Coefficient',' Standard Error','T value','P>|T|','')
  print(object$Coefs)
  cat('(* indicates significnace at the 5% level.)\n\n')
  cat('Observations, Predictions, Residuals, Leverage and Influence\n')
  cat(paste0('Leverage Limit:\t',round(object$LevLim,4),'\n'))
  cat(paste0('Influence Limit:\t',round(object$InflLim,4),'\n'))
  temp <- cbind(object$Y,object$fitted.values,object$ResLevInf)
  temp <- round(temp,4)
  names(temp) <- c('Observation','Prediction','Residual','Leverage','Influence')
  temp$Leverage <- ifelse(object$LevInf.Sig[,1],
                          paste0(temp$Leverage,'*'),temp$Leverage)
  temp$Influence <- ifelse(object$LevInf.Sig[,2],
                           paste0(temp$Influence,'*'),temp$Influence)
  print(temp)
}
#' @rdname print
#' @return \code{print.WREG.WLS} Prints a summary of output list from WREG.WLS
#' @examples 
#' ## print.WREG.WLS
#' # Import some example data
#' peakFQdir <- paste0(
#'   file.path(system.file("exampleDirectory", package = "WREG"),
#'     "pfqImport"))
#' gisFilePath <- file.path(peakFQdir, "pfqSiteInfo.txt")
#' importedData <- importPeakFQ(pfqPath = peakFQdir, gisFile = gisFilePath)
#' 
#' # Organizing input data
#' lp3Data <- importedData$LP3f
#' lp3Data$K <- importedData$LP3k$AEP_0.5
#' Y <- importedData$Y$AEP_0.5
#' X <- importedData$X[c("Sand", "OutletElev", "Slope")]
#' recordLengths <- importedData$recLen
#' transY <- "none"
#' 
#' # Run WLS regression
#' result <- WREG.WLS(Y, X, recordLengths, LP3 = lp3Data, transY)
#' print(result)
#' @export
print.WREG.WLS <- function(x, ...) {
  object <- x
  cat(paste0("Regression Model for ",names(object$Y),'\n'))
  cat('Coefficients fit by weighted least-squares.\n')
  cat("\nCall:\n", object$Inputs$call,'\n\n')
  cat(paste0('\nPerformance Metrics\n(Note: Units are based',
             ' on the transformation of the dependent variable.)\n'))
  cat('Mean Squared Error:\t',
      round(object$PerformanceMetrics$MSE,digits=4),'\n')
  cat(paste0('Root Mean Squared Error:\t',
             round(object$PerformanceMetrics$RMSE,digits=2),'%\n'))
  cat(paste0('Coefficient of Determination (R2):\t',
             round(100*object$PerformanceMetrics$R2,digits=2),'%\n'))
  cat(paste0('Adjusted Coefficient of Determination (Adj-R2):\t',
             round(100*object$PerformanceMetrics$R2_adj,digits=2),'%\n'))
  cat(paste0('Pseudo Coefficient of Determination (Pseudo-R2):\t',
             round(100*object$PerformanceMetrics$R2_pseudo,digits=2),'%\n'))
  cat('Average Variance of Prediction:\t',
      round(object$PerformanceMetrics$AVP,digits=4),'\n')
  cat(paste0('Average Standard Error of Prediction:\t',
             round(object$PerformanceMetrics$Sp,digits=2),'%\n'))
  cat('Model Error Variance:\t',
      round(object$PerformanceMetrics$ModErrVar,digits=4),'\n')
  cat(paste0('Standard Model Error Variance:\t',
             round(object$PerformanceMetrics$StanModErr,digits=2),'%\n\n'))  
  cat('Coefficients of Model')
  object$Coefs[,1:3] <- round(object$Coefs[,1:3],3)
  object$Coefs[,4] <- round(object$Coefs[,4],4)
  object$Coefs[,5] <- ifelse(object$Coefs[,4]<=0.05,'*','')
  names(object$Coefs) <- c('Coefficient',' Standard Error','T value','P>|T|','')
  print(object$Coefs)
  cat('(* indicates significnace at the 5% level.)\n\n')
  cat('Observations, Predictions, Residuals, Leverage and Influence\n')
  cat(paste0('Leverage Limit:\t',round(object$LevLim,4),'\n'))
  cat(paste0('Influence Limit:\t',round(object$InflLim,4),'\n'))
  temp <- cbind(object$Y,object$fitted.values,object$ResLevInf,
                object$PerformanceMetrics$VP.PredVar)
  temp <- round(temp,4)
  names(temp) <- c('Observation','Prediction','Residual','Leverage',
                   'Influence','Variance of Prediction')
  temp$Leverage <- ifelse(object$LevInf.Sig[,1],
                          paste0(temp$Leverage,'*'),temp$Leverage)
  temp$Influence <- ifelse(object$LevInf.Sig[,2],
                           paste0(temp$Influence,'*'),temp$Influence)
  print(temp)
}
#'
#' @rdname print
#' @return \code{print.WREG.GLS} Prints a summary of output list from WREG.GLS
#' @examples
#' ## print.WREG.GLS
#' # Import some example data
#' peakFQdir <- paste0(
#'   file.path(system.file("exampleDirectory", package = "WREG"),
#'     "pfqImport"))
#' gisFilePath <- file.path(peakFQdir, "pfqSiteInfo.txt")
#' importedData <- importPeakFQ(pfqPath = peakFQdir, gisFile = gisFilePath)
#' 
#' # Organizing input data
#' lp3Data <- importedData$LP3f
#' lp3Data$K <- importedData$LP3k$AEP_0.5
#' Y <- importedData$Y$AEP_0.5
#' X <- importedData$X[c("Sand", "OutletElev", "Slope")]
#' recordLengths <- importedData$recLen
#' basinChars <- importedData$BasChars
#' transY <- "none"
#' 
#' # Run GLS regression
#' result <- WREG.GLS(Y, X, recordLengths, LP3 = lp3Data, basinChars, transY)
#' print(result)
#' @export
print.WREG.GLS <- function(x, ...) {
  object <- x
  cat(paste0("Regression Model for ",names(object$Y),'\n'))
  cat('Coefficients fit by generalized least-squares.\n')
  cat("\nCall:\n", object$Inputs$call,'\n\n')
  cat(paste0('\nPerformance Metrics\n(Note: Units are based',
             ' on the transformation of the dependent variable.)\n'))
  cat('Mean Squared Error:\t',
      round(object$PerformanceMetrics$MSE,digits=4),'\n')
  cat(paste0('Root Mean Squared Error:\t',
             round(object$PerformanceMetrics$RMSE,digits=2),'%\n'))
  cat(paste0('Coefficient of Determination (R2):\t',
             round(100*object$PerformanceMetrics$R2,digits=2),'%\n'))
  cat(paste0('Adjusted Coefficient of Determination (Adj-R2):\t',
             round(100*object$PerformanceMetrics$R2_adj,digits=2),'%\n'))
  cat(paste0('Pseudo Coefficient of Determination (Pseudo-R2):\t',
             round(100*object$PerformanceMetrics$R2_pseudo,digits=2),'%\n'))
  cat('Average Variance of Prediction:\t',
      round(object$PerformanceMetrics$AVP,digits=4),'\n')
  cat(paste0('Average Standard Error of Prediction:\t',
             round(object$PerformanceMetrics$Sp,digits=2),'%\n'))
  cat('Model Error Variance:\t',
      round(object$PerformanceMetrics$ModErrVar,digits=4),'\n')
  cat(paste0('Standard Model Error Variance:\t',
             round(object$PerformanceMetrics$StanModErr,digits=2),'%\n\n'))  
  cat('Coefficients of Model')
  object$Coefs[,1:3] <- round(object$Coefs[,1:3],3)
  object$Coefs[,4] <- round(object$Coefs[,4],4)
  object$Coefs[,5] <- ifelse(object$Coefs[,4]<=0.05,'*','')
  names(object$Coefs) <- c('Coefficient',' Standard Error','T value','P>|T|','')
  print(object$Coefs)
  cat('(* indicates significnace at the 5% level.)\n\n')
  cat('Observations, Predictions, Residuals, Leverage and Influence\n')
  cat(paste0('Leverage Limit:\t',round(object$LevLim,4),'\n'))
  cat(paste0('Influence Limit:\t',round(object$InflLim,4),'\n'))
  temp <- cbind(object$Y,object$fitted.values,object$ResLevInf,
                object$PerformanceMetrics$VP.PredVar)
  temp <- round(temp,4)
  names(temp) <- c('Observation','Prediction','Residual','Leverage',
                   'Influence','Variance of Prediction')
  temp$Leverage <- ifelse(object$LevInf.Sig[,1],
                          paste0(temp$Leverage,'*'),temp$Leverage)
  temp$Influence <- ifelse(object$LevInf.Sig[,2],
                           paste0(temp$Influence,'*'),temp$Influence)
  print(temp)
}
#'
#' @rdname print
#' @return \code{print.GLSs} Prints a summary of output list from WREG.GLS with
#' an adjustment for uncertainty in the skewness
#' @examples
#' ## print.WREG.GLS
#' # Import some example data
#' peakFQdir <- paste0(
#'   file.path(system.file("exampleDirectory", package = "WREG"),
#'     "pfqImport"))
#' gisFilePath <- file.path(peakFQdir, "pfqSiteInfo.txt")
#' importedData <- importPeakFQ(pfqPath = peakFQdir, gisFile = gisFilePath)
#' 
#' # Organizing input data
#' lp3Data <- importedData$LP3f
#' lp3Data$K <- importedData$LP3k$AEP_0.5
#' Y <- importedData$Y$AEP_0.5
#' X <- importedData$X[c("Sand", "OutletElev", "Slope")]
#' recordLengths <- importedData$recLen
#' basinChars <- importedData$BasChars
#' transY <- "none"
#' 
#' # Run GLS regression with uncertainty in the skewness
#' result <- WREG.GLS(Y, X, recordLengths, LP3 = lp3Data, basinChars, transY, 
#'   regSkew = TRUE, MSEGR = 0.302)
#' print(result)
#' @export
print.WREG.GLSs <- function(x, ...) {
  object <- x
cat(paste0("Regression Model for ",names(object$Y),'\n'))
cat(paste0('Coefficients fit by generalized least-squares with an \n',
           'adjustment for uncertainty in regional skew.\n'))
cat("\nCall:\n", object$Inputs$call,'\n\n')
cat(paste0('\nPerformance Metrics\n(Note: Units are based',
           ' on the transformation of the dependent variable.)\n'))
cat('Mean Squared Error:\t',
    round(object$PerformanceMetrics$MSE,digits=4),'\n')
cat(paste0('Root Mean Squared Error:\t',
           round(object$PerformanceMetrics$RMSE,digits=2),'%\n'))
cat(paste0('Coefficient of Determination (R2):\t',
           round(100*object$PerformanceMetrics$R2,digits=2),'%\n'))
cat(paste0('Adjusted Coefficient of Determination (Adj-R2):\t',
           round(100*object$PerformanceMetrics$R2_adj,digits=2),'%\n'))
cat(paste0('Pseudo Coefficient of Determination (Pseudo-R2):\t',
           round(100*object$PerformanceMetrics$R2_pseudo,digits=2),'%\n'))
cat('Average Variance of Prediction:\t',
    round(object$PerformanceMetrics$AVP,digits=4),'\n')
cat(paste0('Average Standard Error of Prediction:\t',
           round(object$PerformanceMetrics$Sp,digits=2),'%\n'))
cat('Model Error Variance:\t',
    round(object$PerformanceMetrics$ModErrVar,digits=4),'\n')
cat(paste0('Standard Model Error Variance:\t',
           round(object$PerformanceMetrics$StanModErr,digits=2),'%\n\n'))  
cat('Coefficients of Model')
object$Coefs[,1:3] <- round(object$Coefs[,1:3],3)
object$Coefs[,4] <- round(object$Coefs[,4],4)
object$Coefs[,5] <- ifelse(object$Coefs[,4]<=0.05,'*','')
names(object$Coefs) <- c('Coefficient',' Standard Error','T value','P>|T|','')
print(object$Coefs)
cat('(* indicates significnace at the 5% level.)\n\n')
cat('Observations, Predictions, Residuals, Leverage and Influence\n')
cat(paste0('Leverage Limit:\t',round(object$LevLim,4),'\n'))
cat(paste0('Influence Limit:\t',round(object$InflLim,4),'\n'))
temp <- cbind(object$Y,object$fitted.values,object$ResLevInf,
              object$PerformanceMetrics$VP.PredVar)
temp <- round(temp,4)
names(temp) <- c('Observation','Prediction','Residual','Leverage',
                 'Influence','Variance of Prediction')
temp$Leverage <- ifelse(object$LevInf.Sig[,1],
                        paste0(temp$Leverage,'*'),temp$Leverage)
temp$Influence <- ifelse(object$LevInf.Sig[,2],
                         paste0(temp$Influence,'*'),temp$Influence)
print(temp)
}
