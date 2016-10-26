#' Import Data from Old WREG Format
#' 
#' @description
#' The \code{importWREG} function reads the WREG inputs from a directory 
#' set up for the old WREG program.
#' 
#' @param wregPath A directory that contains all of the files needed to 
#' implement the MatLab version of WREG.
#' @param sites (optional) A vetor of sites that should be return.  Allows for 
#' data subsetting.
#' 
#' @details
#' This functions allows users to use the legacy format of WREG.  This includes
#' an established directory that contains valid \dQuote{SiteInfo.txt}, 
#' \dQuote{FlowChar.txt}, \dQuote{LP3G.txt}, \dQuote{LP3K.txt}, 
#' \dQuote{LP3s.txt} and \dQuote{USGS##########.txt}, multiple files 
#' containing time series for each site.  The file \dQuote{UserWLS.txt} is 
#' optional.  For further information on the format of these files, see the 
#' program manual (Techniques and Methods 4-A8).  Files that are not valid 
#' inputs files for the WREG version describe therein will not be accepted.
#' 
#' @return All outputs are returned as part of a list.  The list includes:
#' \item{sites}{A vector of site IDs.}
#' \item{Y}{A data frame whose comlumns represent unique frequency events, 
#' while the row represent particular sites in the same order as \code{sites}.}
#' \item{X}{A data frame whose columns represent basin characteristics to be used
#'  as dependent variables and whose rows represent sites corresponding 
#'  to \code{sites}.}
#' \item{LP3f}{A matrix containing the fitted LP# parameters that are fixed 
#' across exceedence probability.  These include the standard deviation, 
#' skew and regional skew for each site.}
#' \item{LP3k}{A matrix of the fitted kappa parameters of the LP3 distribution
#'  for each \code{AEP}.}
#' \item{BasChars}{A matrix containing the site IDs, latitudes and longitudes.}
#' \item{recLen}{A square matrix indicating the number of overlapping years
#'  for each site pair.}
#' \item{recCor}{A matrix of the correlaiton between site paris.}
#' \item{UW}{A matrix of user weights, if included.}
#' 
#'@examples
#'wregDir <- file.path(system.file("exampleDirectory", package = "WREG"),
#'  "matlabImport")
#'importedData <- importWREG(wregPath = wregDir)
#'
#'@export
importWREG_General <- function(wregPath,sites='') {
  # Developed by William Farmer, 10 February 2016
  #
  # In this format, the user must indicate the AEP and MSEGR manually.
  
  # wregPath <- file.path('..','SampleInputFiles')
  # sites <- c('5314900','5316900','5316920','5317845','5317850')
  
  # Load and parse SiteInfo.txt
  siteInfoFile <- file.path(wregPath,'SiteInfo.txt')
  if (!file.exists(siteInfoFile)) {
    stop(paste('Could not find',siteInfoFile))
  }
  
  
  # siteInfo <- read.table(siteInfoFile,sep='\t',header=T,
  #                        colClasses = list(Station.ID='character'))
  siteInfo <- read.table(siteInfoFile,sep='\t',header=T,fill=TRUE)
  
  # Check to see if the file has standard column names, otherwise rename the first three columns to the apropriate name
  if (Reduce('&',(is.element(names(siteInfo)[1:25], c('stationID',	'latitude',	'longitude',	'regionalSkew',	'skew',	'standardDeviation',	'flowCharQ2',
                                                'flowCharQ5',	'flowCharQ10',	'flowCharQ25',	'flowCharQ50',	'flowCharQ100',	'flowCharQ200',	'flowCharQ500',
                                                'flowCharQ2.K',	'flowCharQ5.K',	'flowCharQ10.K',	'flowCharQ25.K',	'flowCharQ50.K',	'flowCharQ100.K',	
                                                'flowCharQ200.K',	'flowCharQ500.K',	'zero.1.nonZero.2',	'freqZero',	'Cont.1.PR.2'	)))) == FALSE){
    stop(paste("Please check naming conventions of the general file format"))
  } 
  
  #convert stationID to character class
  siteInfo$stationID <- toupper(as(siteInfo$stationID, 'character'))
  
  #check to see if any sites need a zero appended
  zero_append <- which(nchar(siteInfo$stationID) != 8)
  siteInfo$stationID[zero_append] <- paste0("0",siteInfo$stationID[zero_append])

  BasChars <- siteInfo[,is.element(names(siteInfo),
                                     c('stationID','Lat','Long'))]
  X <- siteInfo[,c(1,26:ncol(siteInfo))]
  sitesOut <- siteInfo$stationID
  Y <- siteInfo[,c(1,7:14)]
  lp3g <- siteInfo[,c(1,5)]
  lp3k <- siteInfo[,c(1,15:23)]
  lp3s <- siteInfo[,c(1,6)]
  
  LP3f <- data.frame(stationID=sitesOut,S=lp3s[,2],G=lp3g[,2],
                     GR=siteInfo$regionalSkew)
  LP3k <- lp3k
  uwls <- NULL
 
  # Load and parse USGS*.txt files
  siteTS <- list()
  
  # See if a USGS Annual time series file exists, if so read it in
  # Otherwise get data from each individual file
  usgs_annual_file = file.path(wregPath,'USGSAnnualTimeSeries.txt')
  if (!file.exists(usgs_annual_file)) {
    stop(paste("Please include a USGS Annual Time Series file with the following name,
               \"USGSAnnualTimeSeries.txt\""))
  }
    
  # split the files in to a list of tables based on site id 
  temp_table <- read.table(usgs_annual_file)
  temp_table[,1] <- toupper(as(temp_table[,1], "character"))
  sites <- unique(temp_table[,1])
  
  for (x in 1:length(sites)){
    siteTS <- c(siteTS, split(temp_table,temp_table[,1] == sites[x])[2])
  }
    
  if (length(siteTS) == 0){
    stop(paste("There are no timeseries files to process.  Check ts naming convetions."))
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
        stop(paste('Overlap in data matrix caused by multiple entries for a year and site.  Revisit data.'))
      }
      
      recCor[i,j] <- recCor[j,i] <- cor(ijdata,jidata)
    }
  }
  
  row.names(recLen) <- row.names(recCor) <- colnames(recLen) <- 
    colnames(recCor) <- sitesOut
  names(LP3k) <- names(Y)
  
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