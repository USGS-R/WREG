#'Weighing Matrix for ROI-WLS (WREG)
#'
#'@description \code{Omega.WLS.ROImatchMatLab} calculates the weighting function
#'for a weighted least-squares regression using regions-of-influence.  This is
#'largely legacy code to match WREG v. 1.05 idiosyncrasies.
#'
#'@param Y.all The dependent variable of interest at all sites in the network,
#'  with any transformations already applied.
#'@param X.all The independent variables  at all sites in the network, with any 
#'  transformations already applied.  Each row represents a site and each column
#'  represents a particular independent variable.  The rows must be in the same
#'  order as the dependent variables in \code{Y}.
#'@param LP3.all A dataframe containing the fitted Log-Pearson Type III standard 
#'  deviate, standard deviation and skew for all sites in the network.  The
#'  names of this data frame are \code{S}, \code{K} and \code{G}.  The order of
#'  the rows must be the same as \code{Y}.
#'@param RecordLengths.all \code{RecordLengths.all} may be a matrix whose rows
#'  and columns are in the same order as \code{Y}.  Each \code{(r,c)} element
#'  represents the length of concurrent record between sites \code{r} and
#'  \code{c}.  The diagonal elements therefore represent each site's full record
#'  length.  For \dQuote{WLS}, only the at-site record lengths are needed.  In
#'  this case, \code{RecordLengths} can be a vector or a matrix.  Here, this
#'  should only include sites in the current region of influence.
#'@param NDX A vector listing the indices of the sites that comprise the region
#'  of influence.
#'  
#'@details This is a legacy function that matches the idiosyncrasies of WREG v.
#'  1.05. This includes using all sites to implement the sigma regression,
#'  averaging across all record lengths, using arbitrary record lengths to
#'  estimate weights and using a new function for estimating the MSE of the
#'  basic OLS model.
#'  
#'  This function will become obsolete once all idiosyncrasies are assessed.
#'  
#'@return \code{Omega.WLS.ROImatchMatLab} returns a list with three elements: 
#'  \item{Omega}{The estimated weighting matrix.} \item{var.modelerror.0}{The
#'  estimated model error variance for a constant-value model.} 
#'  \item{var.modelerror.k}{The estimated model error variance for a k-variable
#'  model.}
#'  
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
#' lp3Data <- importedData$LP3f
#' lp3Data$K <- importedData$LP3k$AEP_0.5
#' Y <- importedData$Y$AEP_0.5
#' X <- importedData$X[c("Sand", "OutletElev", "Slope")]
#' 
#' #### Geographic Region-of-Influence
#' i <- 1 # Site of interest
#' n <- 10 # size of region of influence
#' Gdist <- vector(length=length(Y)) # Empty vector for geographic distances
#' for (j in 1:length(Y)) {
#'   if (i!=j) {
#'     #### Geographic distance
#'     Gdist[j] <- Dist.WREG(Lat1 = importedData$BasChars$Lat[i],
#'       Long1 = importedData$BasChars$Long[i],
#'       Lat2 = importedData$BasChars$Lat[j],
#'       Long2 = importedData$BasChars$Long[j]) # Intersite distance, miles
#'   } else {
#'     Gdist[j] <- Inf # To block self identification.
#'   }
#' }
#' temp <- sort.int(Gdist,index.return=TRUE)
#' NDX <- temp$ix[1:n] # Sites to use in this regression
#' 
#' # Compute weighting matrix
#' weightingResult <- Omega.WLS.ROImatchMatLab(Y.all = Y, X.all = X,
#'   LP3.all = lp3Data, RecordLengths.all = importedData$recLen, NDX = NDX)
#'
#'@export
Omega.WLS.ROImatchMatLab <- function(Y.all,X.all,LP3.all,RecordLengths.all,NDX) {
  # William Farmer, USGS, January 26, 2015
  #
  # WREG v 1.05 (MatLab) is inconsistent in which components are used to develop equations in RoI-WLS
  # Eq. 15 is implemented with all sites
  # The mean inverse record length is implemented with all sites
  
  # Some upfront error handling
  err <- FALSE
  if ((!missing(X.all)&!missing(Y.all))&&
      (length(Y.all)!=nrow(X.all))) {
    warning(paste0("The length of Y.all must be the same as ",
      "the number of rows in X.all."))
    err <- TRUE
  }
  if (missing(Y.all)) {
    warning("Dependent variable (Y.all) must be provided.")
    err <- TRUE
  } else {
    if (!is.numeric(Y.all)) {
      warning("Dependent variable (Y.all) must be provided as class numeric.")
      err <- TRUE
    } else {
      if (sum(is.na(Y.all))>0) {
        warning(paste0("The depedent variable (Y.all) contains missing ",
          "values.  These must be removed."))
        err <- TRUE
      }
      if (sum(is.infinite(Y.all))>0) {
        warning(paste0("The depedent variable (Y.all) contains infinite ",
          "values.  These must be removed."))
        err <- TRUE
      }
    }
  }
  if (missing(X.all)) {
    warning("Independent variables (X.all) must be provided.")
    err <- TRUE
  } else {
    if ((length(unique(apply(X.all,FUN=class,MARGIN=2)))!=1)|
        (unique(apply(X.all,FUN=class,MARGIN=2))!="numeric")) {
      warning("Independent variables (X.all) must be provided as class numeric.")
      err <- TRUE
    } else {
      if (sum(is.na(as.matrix(X.all)))>0) {
        warning(paste0("Some independent variables (X.all) contain missing ",
          "values.  These must be removed."))
        err <- TRUE
      }
      if (sum(is.infinite(as.matrix(X.all)))>0) {
        warning(paste0("Some independent variables (X.all) contain infinite ",
          "values.  These must be removed."))
        err <- TRUE
      }
    }
  }
  if (missing(NDX) | !is.numeric(NDX) | !is.vector(NDX)) {
    warning(paste0("NDX must be provided as a numeric vector."))
    err <- TRUE
  }
  if (!sum(is.element(NDX, 1:length(Y))) == length(NDX)) {
    warning(paste0("NDX must be valid indices of inputs like Y.all"))
    err <- TRUE
  }
  # Error checking LP3
  if (missing(LP3.all)) {
    warning("LP3.all must be provided as an input.")
    err <- TRUE
  } else {
      if (!is.data.frame(LP3.all)) {
        warning(paste("LP3.all must be provided as a data frame with elements named",
          "'S', 'K' and 'G' for standard deivation, deviate and skew,",
          "respectively."))
        err <- TRUE
      } else {
        if (sum(is.element(c("S","K","G"),names(LP3.all)))!=3) {
          warning(paste("In valid elements: The names of the elements in LP3.all are",
            names(LP3.all),". LP3.all must be provided as a data frame with elements named",
            "'S', 'K' and 'G' for standard deivation, deviate and skew,",
            "respectively."))
          err <- TRUE
        }
        if ((length(unique(apply(cbind(LP3.all$S,LP3.all$K,LP3.all$G),FUN=class,MARGIN=2)))!=1)|
            (unique(apply(cbind(LP3.all$S,LP3.all$K,LP3.all$G),FUN=class,MARGIN=2))!="numeric")) {
          warning("LP3.all must be provided as a numeric array")
          err <- TRUE
        } else {
          if (sum(is.infinite(LP3.all$S),is.infinite(LP3.all$K),is.infinite(LP3.all$G))>0) {
            warning(paste0("Some elements of LP3.all$S, LP3.all$K, and LP3.all$G contain infinite ",
              "values.  These must be removed."))
            err <- TRUE
          }
          if (sum(is.na(LP3.all$S),is.na(LP3.all$K),is.na(LP3.all$G))>0) {
            warning(paste0("Some elements of LP3.all$S, LP3.all$K, and LP3.all$G contain missing ",
              "values.  These must be removed."))
            err <- TRUE
          }
        }
      }

  }
  if (missing(RecordLengths.all)) {
    warning("A matrix of RecordLengths.all must be provided as input.")
    err <- TRUE
  } else {
    if (ncol(RecordLengths.all)!=nrow(RecordLengths.all)) {
      warning("RecordLengths.all must be provided as a square array")
      err <- TRUE
    }
    if (!is.numeric(RecordLengths.all)) {
      warning("RecordLengths.all must be provided as a numeric array")
      err <- TRUE
    }
  }
  if (err) {
    stop("Invalid inputs were provided.  See warnings().")
  }
  
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
  B_hat <- solve(t(X)%*%solve(Omega)%*%as.matrix(X))%*%t(X)%*%solve(Omega)%*%Y # OLS estimated coefficients
  Y_hat <- as.matrix(X)%*%B_hat # OLS model estimates
  e <- Y-Y_hat # OLS model residuals
  MSE.OLS <- sum(e^2)/(nrow(X)-ncol(X)) # OLS mean square-error (k-variable)
  MSE.OLS <- sd(e)^2 # BUG: MatLab code uses this expression in RoI WLS.  See line 1450 of WREG v 1.05
  MSE.OLS.0 <- sum((Y-matrix(1,ncol=1,nrow=length(Y))%*%solve(t(matrix(1,ncol=1,nrow=length(Y)))%*%solve(Omega)%*%matrix(1,ncol=1,nrow=length(Y)))%*%t(matrix(1,ncol=1,nrow=length(Y)))%*%solve(Omega)%*%Y)^2)/(nrow(X)-1) # OLS mean square-error (constant model)
  
  ## Estimate k-variable model-error variance
  Omega <- diag(nrow(X.all)) # Temporary weighting matrix (identity) for k-variable sigma regression. Eq 15.
  B.SigReg <- solve(t(X.all)%*%solve(Omega)%*%as.matrix(X.all))%*%t(X.all)%*%solve(Omega)%*%LP3.all$S # OLS estimated coefficients for model of LP3 standard deviations
  # BUG: MatLab fits beta regression (Eq 15) using all sites.  See lines 1305-1316 and 1413-1419 of WREG v 1.05
  Yhat.SigReg <- as.matrix(X.all)%*%B.SigReg # Estimates of sigma regression
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
  Omega <- diag((var.modelerror.k+c1/RecordLengths.all[1:length(NDX)])) # WLS weighting matrix.  Eq 12
  # BUG: MatLab code mis-shuffles the record lengths.  See lines 1457-1459 of WREG v. 1.05
  
  ## Output
  Output <- list(Omega=Omega,var.modelerror.0=var.modelerror.0,var.modelerror.k=var.modelerror.k)
  return(Output)
}
