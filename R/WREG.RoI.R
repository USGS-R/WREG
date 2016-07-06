#'Region-of-Influence Regression (WREG)
#'
#'@description \code{WREG.ROI} implements region-of-influence regression in the
#'WREG framework.
#'
#'@param Y The dependent variable of interest, with any transformations already 
#'  applied.
#'@param X The independent variables in the regression, with any transformations
#'  already applied.  Each row represents a site and each column represents a
#'  particular independe variable.  (If a leading constant is used, it should be
#'  included here as a leading column of ones.)  The rows must be in the same
#'  order as the dependent variables in \code{Y}.
#'@param Reg A string indicating which type of regression should be applied. The
#'  options include: \dQuote{OLS} for ordinary least-squares regression, 
#'  \dQuote{WLS} for weighted least-squares regression, \dQuote{GLS} for
#'  generalized least-squares regression, with no uncertainty in regional skew
#'  and \dQuote{GLSskew} for generalized least-squares regression with
#'  uncertainty in regional skew. (In the case of \dQuote{GLSskew}, the
#'  uncertainty in regional skew must be provided as the mean squared error in
#'  regional skew.)
#'@param transY A required character string indicating if the the 
#'  dependentvariable was transformed by the common logarithm ('log10'), 
#'  transformed by the natural logarithm ('ln') or untransformed ('none').
#'@param recordLengths This input is required for \dQuote{WLS}, \dQuote{GLS} and
#'  \dQuote{GLSskew}.  For \dQuote{GLS} and \dQuote{GLSskew},
#'  \code{recordLengths} should be a matrix whose rows and columns are in the
#'  same order as \code{Y}.  Each \code{(r,c)} element represents the length of
#'  concurrent record between sites \code{r} and \code{c}.  The diagonal
#'  elements therefore represent each site's full record length.  For
#'  \dQuote{WLS}, the only the at-site record lengths are needed. In the case of
#'  \dQuote{WLS}, \code{recordLengths} can be a vector or the matrix described
#'  for \dQuote{GLS} and \dQuote{GLSskew}.
#'@param LP3 A dataframe containing the fitted Log-Pearson Type III standard 
#'  deviate, standard deviation and skew for each site.  The names of this data
#'  frame are \code{S}, \code{K} and \code{G}.  For \dQuote{GLSskew}, the
#'  regional skew value must also be provided in a variable called \code{GR}. 
#'  The order of the rows must be the same as \code{Y}.
#'@param regSkew A logical vector indicating if regional skews are provided with
#'  an adjustment required for uncertainty therein (\code{TRUE}).  The default 
#'  is \code{FALSE}.
#'@param alpha A number, required only for \dQuote{GLS} and \dQuote{GLSskew}. 
#'  \code{alpha} is a parameter used in the estimated cross-correlation between
#'  site records.  See equation 20 in the WREG v. 1.05 manual.  The arbitrary,
#'  default value is 0.01.  The user should fit a different value as needed.
#'@param theta A number, required only for \dQuote{GLS} and \dQuote{GLSskew}. 
#'  \code{theta} is a parameter used in the estimated cross-correlation between
#'  site records.  See equation 20 in the WREG v. 1.05 manual.  The arbitrary,
#'  default value is 0.98.  The user should fit a different value as needed.
#'@param BasinChars A dataframe containing three variables: \code{StationID} is
#'  the numerical identifier (without a leading zero) of each site, \code{Lat}
#'  is the latitude of the site, in decimal degrees, and \code{Long} is the
#'  longitude of the site, in decimal degrees.  The sites must be presented in
#'  the same order as \code{Y}.  Required only for \dQuote{GLS}.
#'@param MSEGR A number. The mean squared error of the regional skew.
#'@param TY A number.  The return period of the event being modeled.  Required
#'  only for \dQuote{GLSskew}.  The default value is \code{2}.  (See the
#'  \code{Legacy} details below.)
#'@param Peak A logical.  Indicates if the event being modeled is a peak flow
#'  event or a low-flow event.  \code{TRUE} indicates a peak flow, while
#'  \code{FALSE} indicates a low-flow event.
#'@param ROI A string indicating how to define the region of influence. 
#'  \dQuote{PRoI} signifies physiographic, independent or predictor-variable
#'  region of influence. \dQuote{GRoI} calls for a geographic region of
#'  influence. \dQuote{HRoI} requests a hybrid region of influence.  Details on
#'  each approach are provided in the manual for WREG v. 1.0.
#'@param n The number of sites to include in the region of influence.
#'@param D Required for \dQuote{HRoI}, the geographic limit within which to
#'  search for a physiographic region of influence.  In WREG v. 1.05 (see
#'  \code{Legacy} below), this is interpretted in meters.  Elsewise, this is
#'  interpretted as miles.
#'@param DistMeth Required for \dQuote{GLS} and \dQuote{GLSskew}.  A value of
#'  \code{1} indicates that the "Nautical Mile" approximation should be used to
#'  calculate inter-site distances.  A value of \code{2} designates the
#'  Haversine approximation.  See \code{\link{Dist.WREG}}.  The default value is
#'  \code{2}.  (See the \code{Legacy} details below.)
#'@param Legacy A logical.  A value of \code{TRUE} forces the WREG program to
#'  behave identically to WREG v. 1.05, with BUGS and all.  It will force
#'  \code{TY=2} and \code{DistMeth=1}.  For ROI regressions, it will also
#'  require a specific calculation for weighing matrices in \dQuote{WLS}
#'  (\code{\link{Omega.WLS.ROImatchMatLab}}), \dQuote{GLS}, and \dQuote{GLSskew}
#'  (see \code{\link{Omega.GLS.ROImatchMatLab}}). \code{Legacy} also forces the
#'  distance limit \code{D} to be interpretted in meters.
#'  
#'@details The support for region-of-influence regression is described in the
#'  manual of WREG v. 1.0.  \code{WREG.RoI} iterates through the sites of
#'  \code{Y}, defines a region of influence and implements the specified
#'  regression by calling a WREG function.
#'  
#'  The logical handle \code{Legacy} has been included to test that this program
#'  returns the same results as WREG v. 1.05.  In the development of this code,
#'  some idiosyncrasies of the MatLab code for WREG v. 1.05 became apparent. 
#'  Setting \code{Legacy} equal to \code{TRUE} forces the code to use the same
#'  idiosycrasies as WREG v. 1.05.  Some of these idosyncrasies may be bugs in
#'  the code.  Further analysis is needed.  For information on the specific
#'  idiosyncrasies, see the notes for the \code{Legacy} input and the links to
#'  other functions in this package.
#'  
#'@return As with other WREG functions, \code{WREG.RoI} returns a large list
#'  of regression parameters and metrics.  This list varies depending on the
#'  \code{Reg} specified, but may contain: \item{fitted.values}{A vector of
#'  model estimates from the regression model.} \item{residuals}{A vector of
#'  model residuals.} \item{PerformanceMetrics}{A list of four elements.  These
#'  represent approximate performance regression across all of the
#'  region-of-influence regressions. These include the mean squared error of
#'  residuals (\code{MSE}), the coefficient of determination (\code{R2}), the
#'  adjusted coefficient of determination (\code{R2_adj}) and the root mean
#'  squared error (\code{RMSE}, in percent).  Details on the appropriateness and
#'  applicability of performance metrics can be found in the WREG manual.} 
#'  \item{Coefficients}{A list composed of four elements: (1) \code{Values}
#'  contains the regression coefficeints estimated for the model built around
#'  each observation, (2) \code{StanError} contains the standard errors
#'  of each regression coefficient for the ROI regressions around each
#'  observations, (3) \code{TStatistic} contains the Student's T-statistic of
#'  each regression coefficient for the ROI regression built around each
#'  observation and (4) \code{pValue} contains the significance probability of
#'  each regression coefficient for the ROI regressions built around each 
#'  observation.  Each element of the list is a matrix the same size as
#'  \code{X}} \item{ROI.Regressions}{A list of elements and outputs from each
#'  individual ROI regression.  These include a matrix of the sites used in each
#'  ROI regression (\code{Sites.Used}, a \code{length(Y)}-by-\code{n} matrix), a
#'  matrix of the geographic distances between the selected sites in each ROI
#'  regression (\code{ Gdist.Used}, a \code{length(Y)}-by-\code{n} matrix), a
#'  matrix of the physiographic distances between the selected sites in each ROI
#'  regression (\code{ Pdist.Used}, a \code{length(Y)}-by-\code{n} matrix), a
#'  matrix of the observations used in each ROI regression (\code{Obs.Used}, a
#'  \code{length(Y)}-by- \code{n} matrix), a matrix of model fits in the region
#'  of influence (\code{Fits}, a \code{length(Y)}-by-\code{n} matrix), a matrix
#'  of model residuals in the region of influence (\code{Residuals}, a
#'  \code{length(Y)}-by-\code{n} matrix), a matrix of leverages for each ROI
#'  regression (\code{Leverage}, a \code{length(Y)}-by-\code{n} matrix), a
#'  matrix of influences for each ROI regression (\code{Influence}, a \code{ 
#'  length(Y)}-by-\code{n} matrix), a logical matrix indicating if the leverage
#'  is significant in each ROI regression (\code{Leverage.Significance}, a
#'  \code{length(Y)} -by-\code{n} matrix), a logical matrix indicating if the
#'  influence is significant in each ROI regression
#'  (\code{Influence.Significance}, a \code{length(Y)} -by-\code{n} matrix), a
#'  vector of critical leverage values for each ROI regression 
#'  (\code{Leverage.Limits}), a vector of critical influence values for each ROI
#'  regression (\code{Influence.Limits}) and a list of performance metrics for
#'  each ROI regression (\code{PerformanceMetrics}).  The last element, 
#'  \code{PerformanceMetrics} is identical to the same output from other 
#'  functions excpet that every element is multiplied by the number of
#'  observations so as to capture the individual performance of each ROI
#'  regression.} \item{ROI.InputParameters}{A list of input parameters to record
#'  the controls on the ROI regression.  \code{D} idicates the limit used in
#'  \dQuote{HRoI}. \code{n} indicates the size of the region of influence. 
#'  \code{ROI} is a string indicating the type of region of influence. 
#'  \code{Legacy} is a logical indicating if the WREG v. 1.05 idiosycrasies were
#'  implemented.}
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
#' X <- importedData$X[c("A")]
#' transY <- "none"
#' basinChars <- importedData$BasChars
#' #result <- WREG.OLS(Y, X, transY)
#' 
#' result <- WREG.RoI(Y, X, Reg = "OLS", transY, BasinChars = basinChars,
#'   ROI='GRoI', n = 10L)
#'@export
WREG.RoI <- function(Y,X,Reg,transY=NA,
  recordLengths = NA,LP3 = NA,regSkew=FALSE,
  alpha=0.01,theta=0.98,BasinChars=NA,MSEGR=NA,TY=2,Peak=T,
  ROI=c('PRoI','GRoI','HRoI'),n=NA,D=250,DistMeth=2,Legacy=FALSE) {
  # William Farmer, USGS, January 23, 2015
  
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
  if (missing(n)|!is.integer(n)|length(n)>1) {
    warning(paste0("n must be provided as a single integer value."))
    err <- TRUE
  }
  # Error checking LP3
  if (is.element(Reg,c("GLS","GLSskew","WLS"))) {
    if (missing(LP3)) {
      warning("LP3 must be provided as an input.")
      err <- TRUE
    } else {
      if (!regSkew) {
        if (!is.data.frame(LP3)) {
          warning(paste("LP3 must be provided as a data frame with elements named",
            "'S', 'K' and 'G' for standard deivation, deviate and skew,",
            "respectively."))
          err <- TRUE
        } else {
          if (sum(is.element(c("S","K","G"),names(LP3)))!=3) {
            warning(paste("In valid elements: The names of the elements in LP3 are",
              names(LP3),". LP3 must be provided as a data frame with elements named",
              "'S', 'K' and 'G' for standard deivation, deviate and skew,",
              "respectively."))
            err <- TRUE
          }
          if ((length(unique(apply(cbind(LP3$S,LP3$K,LP3$G),FUN=class,MARGIN=2)))!=1)|
              (unique(apply(cbind(LP3$S,LP3$K,LP3$G),FUN=class,MARGIN=2))!="numeric")) {
            warning("LP3 must be provided as a numeric array")
            err <- TRUE
          } else {
            if (sum(is.infinite(LP3$S),is.infinite(LP3$K),is.infinite(LP3$G))>0) {
              warning(paste0("Some elements of LP3$S, LP3$K, and LP3$G contain infinite ",
                "values.  These must be removed."))
              err <- TRUE
            }
            if (sum(is.na(LP3$S),is.na(LP3$K),is.na(LP3$G))>0) {
              warning(paste0("Some elements of LP3$S, LP3$K, and LP3$G contain missing ",
                "values.  These must be removed."))
              err <- TRUE
            }
          }
        }
      } else {
        if (!is.data.frame(LP3)) {
          warning(paste("LP3 must be provided as a data frame with elements named",
            "'S', 'K', 'G' and 'GR' for standard deivation, deviate,",
            "skew and regional skew, respectively."))
          err <- TRUE
        } else {
          if (sum(is.element(c("S","K","G","GR"),names(LP3)))!=4) {
            warning(paste("In valid elements: The names of the elements in LP3 are",
              names(LP3),". LP3 must be provided as a data frame with elements named",
              "'S', 'K', 'G' and 'GR' for standard deivation, deviate,",
              "skew and regional skew, respectively."))
            err <- TRUE
          }
          if ((length(unique(apply(cbind(LP3$S,LP3$K,LP3$G,LP3$GR),FUN=class,MARGIN=2)))!=1)|
              (unique(apply(cbind(LP3$S,LP3$K,LP3$G,LP3$GR),FUN=class,MARGIN=2))!="numeric")) {
            warning("LP3 must be provided as a numeric array")
            err <- TRUE
          } else {
            if (sum(is.infinite(LP3$S),is.infinite(LP3$K),
              is.infinite(LP3$G),is.infinite(LP3$GR))>0) {
              warning(paste0("Some elements of LP3$S, LP3$K, LP3$G and LP3$GR contain ",
                "infinite values.  These must be removed."))
              err <- TRUE
            }
            if (sum(is.na(LP3$S),is.na(LP3$K),is.na(LP3$G),is.na(LP3$GR))>0) {
              warning(paste0("Some elements of LP3$S, LP3$K, LP3$G and LP3$GR contain ",
                "missing values.  These must be removed."))
              err <- TRUE
            }
          }
        }
      }
    }
    if (missing(recordLengths)) {
      warning("A matrix of recordLengths must be provided as input.")
      err <- TRUE
    } else {
      if (ncol(recordLengths)!=nrow(recordLengths)) {
        warning("recordLengths must be provided as a square array")
        err <- TRUE
      }
      if (!is.numeric(recordLengths)) {
        warning("recordLengths must be provided as a numeric array")
        err <- TRUE
      }
    }
  }
  if (missing(BasinChars)) {
    warning("BasinChars must be provided as input.")
    err <- TRUE
  } else {
    if (!is.data.frame(BasinChars)) {
      warning(paste("'BasinChars' must be provided as a data frame with elements",
        "named 'Station.ID', 'Lat' and 'Long' for standard deivation,",
        "deviate and skew, respectively."))
      err <- TRUE
    } else {
      if (sum(is.element(c("Station.ID","Lat","Long"),names(BasinChars)))!=3) {
        warning(paste("In valid elements: The names of the elements in",
          "BasinChars are",names(BasinChars),
          ".  'BasinChars' must be provided as a data frame with elements",
          "named 'Station.ID', 'Lat' and 'Long' for standard deivation,",
          "deviate and skew, respectively."))
        err <- TRUE
      } else {
        if ((length(unique(apply(cbind(BasinChars$Lat,BasinChars$Long),FUN=class,MARGIN=2)))!=1)|
            (unique(apply(cbind(BasinChars$Lat,BasinChars$Long),FUN=class,MARGIN=2))!="numeric")) {
          warning("latitudes and longitudes must be provided as class numeric.")
          err <- TRUE
        } else {
          if (sum(is.na(c(BasinChars$Lat,BasinChars$Long)))>0) {
            warning(paste0("Some latitudes and longitudes contain missing ",
              "values.  These must be removed."))
            err <- TRUE
          }
          if (sum(is.infinite(c(BasinChars$Lat,BasinChars$Long)))>0) {
            warning(paste0("Some latitudes and longitudes contain infinite ",
              "values.  These must be removed."))
            err <- TRUE
          }
        }
      }
    }
  }
  
  ## Add controls to meet legacy demands
  if (!is.logical(Legacy)) {
    warning("Legacy must be either TRUE to force matching with previous",
      "versions or FALSE for correct computations.")
    err <- TRUE
  }
  ##    NOTE: Legacy forces program to return the same results as WREG v 1.05.
  if (Legacy) { # If legacy is indicated, override custom inputs.
    TY <- 2 # WREG v1.05 does not read this input correctly.
    DistMeth <- 1 # WREG v1.05 uses "Nautical Mile" approximation
  }
  
  if (err) {
    stop("Invalid inputs were provided.  See warnings().")
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
    x0.i.nocon <- t(matrix(unlist(rep(X.nocon[i,],times=length(Y))),ncol=length(Y))) # Predictors at target site, less a constant
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
    if (!is.na(sum(c(recordLengths)))&&is.matrix(recordLengths)) {# GLS
      recordLengths.i <- recordLengths[NDX,NDX] # Record lengths from region of influence
    } else if (!is.na(sum(c(recordLengths)))&&is.vector(recordLengths)) { # WLS
      recordLengths.i <- recordLengths[NDX] # Record lengths from region of influence
    } else {
      recordLengths.i <- recordLengths
    }
    BasinChars.i <- BasinChars[NDX,] # Basin characteristics (IDs, Lat, Long) from region of influence.
    LP3.i <- data.frame(LP3)[NDX,] # LP3 parameters from region of influence
    if (Legacy) { # Apply regressions to match MatLab WREG v 1.05.
      if (Reg=='WLS') { # Use subroutine to correct WLS weights to meet v1.05 idiosyncrasies
        WeightFix <- Omega.WLS.ROImatchMatLab(Y.all=Y,X.all=X,LP3.all=LP3,RecordLengths.all=recordLengths,NDX=NDX)
        Reg.i <- WREG.UW(Y = Y.i, 
          X = as.matrix(X.i , ncol = length(X.i) / length(Y)),
          customWeight = WeightFix, transY, x0 = x0.i)
      } else if (is.element(Reg,c('GLS','GLSskew'))) { # Use subroutine to correct GLS and GLSskew weights to meet v1.05 idiosyncrasies
        # "Corrected" weighting matrix to match MATLAB code.  Returns weighting matrix and var.moderror.k
        WeightFix.k <- Omega.GLS.ROImatchMatLab(alpha=alpha,theta=theta,Independent=BasinChars.i,X=X.i,Y=Y.i,RecordLengths=recordLengths.i,LP3=LP3.i,MSEGR=MSEGR,TY=TY,Peak=Peak,X.all=X,LP3.all=LP3,DistMeth=DistMeth)
        # "Corrected" to match MATLAB code.  Returns weighting matrix and var.moderror.0
        WeightFix.0 <- Omega.GLS.ROImatchMatLab(alpha=alpha,theta=theta,Independent=BasinChars.i,X=matrix(1,ncol=1,nrow=nrow(X.i)),Y=Y.i,RecordLengths=recordLengths.i,LP3=LP3.i,MSEGR=MSEGR,TY=TY,Peak=Peak,X.all=matrix(1,ncol=1,nrow=nrow(X)),LP3.all=LP3,DistMeth=DistMeth)
        WeightFix <- list(Omega=WeightFix.k$Omega,var.modelerror.0=WeightFix.0$GSQ,var.modelerror.k=WeightFix.k$GSQ)
        Reg.i <- WREG.UW(Y = Y.i, 
          X = as.matrix(X.i , ncol = length(X.i) / length(Y)),
          customWeight = WeightFix, transY, x0 = x0.i)
      } else if (Reg=='OLS') { # Apply OLS normally.
        Reg.i <- WREG.OLS(Y = Y.i, 
          X = as.matrix(X.i , ncol = length(X.i) / length(Y)), 
          transY, x0 = x0.i)
      }
    } else { # Aply WREG with "bugs" corrected.
      if (Reg=='WLS') {
        Reg.i <- WREG.WLS(Y = Y.i, 
          X = as.matrix(X.i , ncol = length(X.i) / length(Y)),
          recordLengths = recordLengths.i,
          LP3 = LP3.i, transY, x0 = x0.i)
      } else if (is.element(Reg,c('GLS','GLSskew'))) { 
        Reg.i <- WREG.GLS(Y = Y.i, 
          X = as.matrix(X.i , ncol = length(X.i) / length(Y)),
          recordLengths = recordLengths.i,
          LP3 = LP3.i, transY,
          x0=x0.i, alpha, theta, Peak, distMeth = DistMeth,
          regSkew, MSEGR, TY, legacy=Legacy)
      } else if (Reg=='OLS') {
        Reg.i <- WREG.OLS(Y = Y.i, 
          X = as.matrix(X.i , ncol = length(X.i) / length(Y)), 
          transY, x0 = x0.i)
      }
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