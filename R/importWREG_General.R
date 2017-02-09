#'Import Data from Generic (v1.06) Files
#'
#'@description The \code{importWREG_General} function reads the WREG inputs from
#'a directory set up with new generic file formats.
#'
#'@param wregPath A directory that contains all of the files needed to implement
#'  the MatLab version of WREG.
#'  
#'@details This function allows users to use a more streamlined data format: 
#'  only two or three files are required.  This includes the 
#'  \dQuote{SiteInfo.txt}, \dQuote{USGSAnnualTimeSeries.txt}, and, optionally, 
#'  \dQuote{UserWLS.txt}.  The \dQuote{SiteInfo.txt} file should contain the
#'  following columns (with headers in parentheses): the station identification
#'  number (stationID), the latitude and longitude of the stations (latitude and
#'  longitude), the regional skew computed for each site (regionalSkew), the 
#'  at-site skew value (skew), the standard deviation of the Log-Pearson 
#'  Type-III distribution used to fit the series (standardDeviation), a series 
#'  of streamflow characteristics to be evaluated (Q#, where # indicates a 
#'  specific return period), Q#.k (the characteristic-specific kappa value from 
#'  the fitted Log-Pearson Type-III distribution, where again # specifies a 
#'  certain return period), and any explanatory variables to be used for 
#'  analysis.  The \dQuote{USGSAnnualTimeSeries.txt} and \dQuote{UserWLS.txt}
#'  files follow the format outlined in USGS Techniques and Methods 4-A8.
#'  
#'@return All outputs are returned as part of a list.  The list includes: 
#'  \item{sites}{A vector of site IDs.} \item{Y}{A data frame whose comlumns
#'  represent unique frequency events, while the row represent particular sites
#'  in the same order as \code{sites}.} \item{X}{A data frame whose columns
#'  represent basin characteristics to be used as dependent variables and whose
#'  rows represent sites corresponding to \code{sites}.} \item{LP3f}{A matrix
#'  containing the fitted LP# parameters that are fixed across exceedence
#'  probability.  These include the standard deviation, skew and regional skew
#'  for each site.} \item{LP3k}{A matrix of the fitted kappa parameters of the
#'  LP3 distribution for each \code{AEP}.} \item{BasChars}{A matrix containing
#'  the site IDs, latitudes and longitudes.} \item{recLen}{A square matrix
#'  indicating the number of overlapping years for each site pair.} 
#'  \item{recCor}{A matrix of the correlaiton between site paris.} \item{UW}{A
#'  matrix of user weights, if included.}
#'  
#'@examples
#'wregDir <- file.path(system.file("exampleDirectory", package = "WREG"),
#'  "generalImport")
#'importedData <- importWREG_General(wregPath = wregDir)
#'
#'@export
importWREG_General <- function(wregPath) {
  # Developed by William Farmer, 10 February 2016
  #
  # In this format, the user must indicate the AEP and MSEGR manually.
  
  # wregPath <- file.path('..','SampleInputFiles')
  
  # Load and parse SiteInfo.txt
  siteInfoFile <- file.path(wregPath,'SiteInfo.txt')
  if (!file.exists(siteInfoFile)) {
    stop(paste('Could not find',siteInfoFile))
  }
  
  siteInfo <- read.table(siteInfoFile,sep='\t',header=T,fill=TRUE)
  
  # Check to see if the file has standard column names, 
  # otherwise rename the first three columns to the apropriate name
  if (sum(is.element(names(siteInfo), c('stationID',	'latitude',	'longitude',
    'regionalSkew',	'skew',	'standardDeviation'))) != 6) {
    stop("Please check naming conventions of the general file format")
  }
  
  #convert stationID to character class
  siteInfo$stationID <- toupper(as(siteInfo$stationID, 'character'))
  
  #check to see if any sites need a zero appended
  zero_append <- which(nchar(siteInfo$stationID) != 8)
  siteInfo$stationID[zero_append] <- paste0("0",siteInfo$stationID[zero_append])

  # Grab basic basin characteristics
  BasChars <- siteInfo[, c('stationID',	'latitude',	'longitude')]
  names(BasChars) <- c('Station.ID', 'Lat', 'Long')
  sitesOut <- siteInfo$stationID
  
  # Check if flow characteristics have been provided correctly
  flowChars <- grep(pattern = "Q[[:digit:]]+", x = names(siteInfo)) != 
    grep(pattern = "Q[[:digit:]]+.k", x = names(siteInfo))
  flowChars <- grep(pattern = "Q[[:digit:]]+", x = names(siteInfo), 
    value = TRUE)[flowChars]
  if (length(flowChars) == 0) {
    stop(paste("Either no flow characteristics were provided or the table", 
      "headings in SiteInfo.txt were incorrect.  Please label flow", 
      "characteristics as Q# where # represents a return period."))
  }
  Y <- siteInfo[, c("stationID", flowChars)]
  
  # Check if kappa values were provided, make LP3 files
  kappas <- grep(pattern = "Q[[:digit:]]+.k", x = names(siteInfo), 
    value = TRUE)
  if (length(kappas) == 0) {
    stop(paste("Either no kappa values were provided or the table headings",
      "in SiteInfo.txt were incorrect.  Please label kappa values as Q#.k",
      "where # represents a return period."))
  }
  LP3k <- siteInfo[, c("stationID", kappas)]
  names(LP3k)[1] <- "Station.ID"
  kappas <- grep(pattern = "Q[[:digit:]]+.k", x = names(LP3k))
  names(LP3k)[kappas] <- unlist(lapply(names(LP3k)[kappas],
    FUN = strsplit, split = ".k"))
  LP3f <- siteInfo[, c("stationID", "standardDeviation", "skew",
    "regionalSkew")]
  names(LP3f) <- c("Station.ID", "S", "G", "GR")
  
  # Check for explanatory variables
  X <- which(!is.element(names(siteInfo), c('stationID',	'latitude',
    'longitude',	'regionalSkew',	'skew',	'standardDeviation', 
    grep(pattern = "Q[[:digit:]]+", x = names(siteInfo), value = TRUE))))
  if (length(X) == 0) {
    stop("No explanatory variables were found in the SiteInfo.txt file.")
  }
  X <- siteInfo[,c("stationID", names(siteInfo)[X])]
  names(X)[1] <- 'Station.ID'
  
  # Load and parse UserWLS.txt (if available)
  uwlsFile <- file.path(wregPath,'UserWLS.txt')
  uwls <- NULL
  if (file.exists(uwlsFile)) {
    uwls <- read.table(uwlsFile,sep='\t',header=F)
    if (nrow(uwls)!=ncol(uwls)) {
      stop('UserWLS.txt must be a square matrix.')
    }
    if (!(identical(site1,site2)&identical(site1,site3a)&
        identical(site1,site3b)&identical(site1,site3c))) {
      uwls <- matrix(NA,ncol=length(site1),nrow=length(site2))
      warning(paste0('Because the input files have a different order of sites,',
        ' the user-provided weighting matrix is invalid.',
        ' It has been replaced with NAs.'))
    } else if (nrow(uwls)!=n2) {
      uwls <- uwls[ndx,ndx]
      warning(paste0('UserWLS.txt was ammended to match sites in arguments.',
        ' (Note: Be sure weights do not depend on sites that were removed.)'))
    }
    uwls <- as.matrix(uwls)
    row.names(uwls) <- colnames(uwls) <- sitesOut
  }
 
  # Load and parse time series files
  siteTS <- list()
  # See if a USGS Annual time series file exists, if so read it in
  # Otherwise get data from each individual file
  usgs_annual_file = file.path(wregPath,'USGSAnnualTimeSeries.txt')
  if (!file.exists(usgs_annual_file)) {
    stop(paste("Please include a USGS Annual Time Series file with",
      "the following name: \"USGSAnnualTimeSeries.txt\""))
  }
    
  # split the files in to a list of tables based on site id 
  temp_table <- read.table(usgs_annual_file)
  temp_table[,1] <- toupper(as(temp_table[,1], "character"))
  sites <- unique(temp_table[,1])
  
  for (x in 1:length(sites)){
    siteTS <- c(siteTS, split(temp_table,temp_table[,1] == sites[x])[2])
  }
    
  if (length(siteTS) == 0){
    stop(paste("There are no timeseries files to process.  Check ts naming", 
      "convetions."))
  }
  
  recLen <- recCor <- matrix(NA,ncol=length(siteTS),nrow=length(siteTS))
  
  for (i in 1:length(siteTS)) {
    recLen[i,i] <- nrow(siteTS[[i]])
    idata <- abs(siteTS[[i]][,2:3])
    for (j in 1:i) {
      jdata <- abs(siteTS[[j]][,2:3])
      overlap <- intersect(idata[,1],jdata[,1])
      recLen[i,j] <- recLen[j,i] <- length(overlap)
      if (length(overlap)==0) {next}
      ijdata <- idata[which(is.element(idata[,1],overlap)),2]
      jidata <- jdata[which(is.element(jdata[,1],overlap)),2]
      if (length(unique(ijdata))==1|length(unique(jidata))==1) {next}
      
      #check to see if the data has the proper dimensions to be correlated
      if (length(ijdata) != length(jidata)){
        stop(paste('Overlap in data matrix caused by multiple entries',
          'for a year and site.  Revisit data.'))
      }
      
      recCor[i,j] <- recCor[j,i] <- cor(ijdata,jidata)
    }
  }
  
  row.names(recLen) <- row.names(recCor) <- colnames(recLen) <- 
    colnames(recCor) <- sitesOut
  
  # Output result
  result <- list(
    sites=sitesOut,
    Y=Y,
    X=X,
    LP3f=LP3f,
    LP3k=LP3k,
    BasChars=BasChars,
    recLen=recLen,
    recCor=recCor,
    UW=uwls
  )
  return(result)
  
}