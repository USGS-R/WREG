#' Import Data from PeakFQ Output
#' 
#' @description
#' The \code{importPeakFQ} function reads output directly from PeakFQ and 
#' prepares it for use in WREG-R.
#' 
#' @param pfqPath A directory that contains all of PeakFQ files for each site 
#' in the GIS file.  Files can be contained in the main path or in 
#' subdirectories of \code{pfqPath}.  Each site should be represented by one
#'  and only one .EXP and .PRT file.
#' @param gisFile A tab-delimited text file that contains a matrix whose rows 
#' represent the sites and contains columns that include \sQuote{Station.ID}, 
#' \sQuote{Lat}, \sQuote{Long} and any other variables to be used for analysis.
#' @param sites (optional) A vetor of sites that should be return.  Allows for 
#' data subsetting.
#' 
#' @details
#' This functions allows users to read output directly from PeakFQ without 
#' the need to manipulate additional files.  The site selection is driven 
#' by the sites in the GIS table.
#' 
#' @return All outputs are returned as part of a list.  The list includes:
#' \item{sites}{A vector of site IDs.}
#' \item{Y}{A matrix whose comlumns represent unique frequency events, 
#' while the row represent particular sites in the same order as \code{sites}.}
#' \item{AEP}{A vector of annual exceedance probabilities associated with the 
#' columns of \code{Y}.}
#' \item{X}{A matrix whose columns represent basin characteristics to be used
#'  as dependent variables and whose rows represent sites corresponding to \code{sites}.}
#' \item{LP3f}{A matrix containing the fitted LP# parameters that are fixed 
#' across exceedence probability.  These include the standard deviation, 
#' skew and regional skew for each site.}
#' \item{LP3k}{A matrix of the fitted kappa parameters of the LP3 distribution
#'  for each \code{AEP}.}
#' \item{BasChars}{A matrix containing the site IDs, latitudes and longitudes.}
#' \item{MSEGR}{The mean-squared-error in the regional skew.}
#' \item{recLen}{A square matrix indicating the number of overlapping years
#'  for each site pair.}
#' \item{recCor}{A matrix of the correlaiton between site paris.}
#' 
#' @examples
#' peakFQdir <- paste0(
#'   file.path(system.file("exampleDirectory", package = "WREG"),
#'     "pfqImport"))
#' gisFilePath <- file.path(peakFQdir, "pfqSiteInfo.txt")
#' importedData <- importPeakFQ(pfqPath = peakFQdir, gisFile = gisFilePath)
#' 
#'@export
importPeakFQ <- function(pfqPath,gisFile,sites='') {
  # Developed by William Farmer, 04 February 2016
  
  # Load GIS file
  gisData <- read.table(file=gisFile,sep='\t',header=T)
  
  # Check to see if the file has standard column names, otherwise rename the first three columns to the apropriate name
  if (Reduce('&',(is.element(names(gisData), c('Station.ID', 'Lat', 'Long')))) == FALSE){
    names(gisData)[1:3] <- c('Station.ID', 'Lat', 'Long')
  } 
  
  #convert Station.ID to character class
  gisData$Station.ID <- as(gisData$Station.ID, 'character')
  
  gisData$Station.ID <- ifelse(nchar(gisData$Station.ID)%%2>0,
    paste0("0",gisData$Station.ID), gisData$Station.ID)
  
  # Remove old WREG variables that are reproduced elsewhere
  remNDX <- -which(is.element(names(gisData), c('No..Annual.Series',
    'Zero.1.NonZero.2',
    'FreqZero',
    'Regional.Skew',
    'Cont.1.PR.2')))
  gisData <- gisData[, remNDX]
  
  # subset to specified sites
  if (sites!='') {
    gisData <- gisData[gisData$Station.ID %in% sites,]
  }
  #ndx <- which(is.element(gisData$Station.ID,sites))
  
  # Independent variables
  BasChars <- gisData[c("Station.ID","Lat","Long")]
  #X <- gisData[ndx,!is.element(names(gisData),
  #                             c('Lat','Long'))]
  
  # Search for EXP for each site
  allFilesEXP <- list.files("*\\.EXP$",
                            path=pfqPath,recursive=T)
  
  ##Look for sites listed in the GIS file in the allFilesEXP vector
  gisSites <- grep(paste(gisData$Station.ID,collapse="|"), 
                   allFilesEXP,
                   ignore.case=TRUE
  )
  
  ##remove files for sites that are not in GIS file
  allFilesEXP <- allFilesEXP[gisSites]
  
  # NEED  a method to handle duplicates
  allFilesEXP <- file.path(pfqPath,allFilesEXP)
  
  # Pulls from EXP curtesy of Janet Curran
  ## Could improve by looking for flag rather than a hard-coded skip
  
  
  EXP_SiteID <- do.call(rbind,lapply(allFilesEXP,read.table,skip=1,nrows=1,colClasses="character"))
  
  ###Send warning message for dropped sites
  ###Get sites that are not in GIS file for warning message
  droppedSites <- gisData$Station.ID[which(!(gisData$Station.ID %in% EXP_SiteID[,2]))]
  
  if(length(droppedSites > 0))
  {
    warning(paste("The following sites are present in the GIS file and not found in EXP files in the Peak FQ output directory:",
                  droppedSites,sep=",")
    )
  }
  
  EXP_G <- do.call(rbind,lapply(allFilesEXP,read.table,skip=7,nrows=1))
  EXP_S <- do.call(rbind,lapply(allFilesEXP,read.table,skip=9,nrows=1))
  EXP_GR <- do.call(rbind,lapply(allFilesEXP,read.table,skip=14,nrows=1))
  EXP_MSEGR <- do.call(rbind,lapply(allFilesEXP,read.table,skip=15,nrows=1))
  EXP_N <- do.call(rbind,lapply(allFilesEXP,read.table,skip=16,nrows=1))[,2] + 
    do.call(rbind,lapply(allFilesEXP,read.table,skip=17,nrows=1))[,2]
  EXP_AEP <- do.call(rbind,lapply(allFilesEXP,read.table,skip=21,nrows=1))
  EXP_Est <- do.call(rbind,lapply(allFilesEXP,read.table,skip=22,nrows=1))
  EXP_Var <- do.call(rbind,lapply(allFilesEXP,read.table,skip=23,nrows=1))
  EXP_K <- do.call(rbind,lapply(allFilesEXP,read.table,skip=26,nrows=1))
  
  # Dependent Variables (with names)
  Y <- EXP_Est[,2:ncol(EXP_Est)]
  AEP <- EXP_AEP[,2:ncol(EXP_AEP)]
  # LP3 Variables
  LP3f <- data.frame(S=EXP_S[,2],G=EXP_G[,2],GR=EXP_GR[,2])
  LP3k <- EXP_K[,2:ncol(EXP_K)]
  
  # Search for PRT for each site
  allFilesPRT <- list.files("*\\.PRT$",
                            path=pfqPath,recursive=T)
  
  ##Look for sites listed in the GIS file in the allFilesPRT vector
  gisSites <- grep(paste(gisData$Station.ID,collapse="|"), 
                   allFilesPRT,
                   ignore.case=TRUE
  )
  
  
  
  ##remove files for sites that are not in GIS file
  allFilesPRT <- allFilesPRT[gisSites]
  
  # NEED  a method to handle duplicates
  allFilesPRT <- file.path(pfqPath,allFilesPRT)
  
  #Function to format the parsed timeseries
  tsFormat <- function(x) {
    x <- x[[2]]
    out <- data.frame(matrix(NA,ncol=2,nrow=length(x)))
    for (i in 1:length(x)) {
      iTemp <- read.table(text=x[[i]],stringsAsFactors=F)
      while (length(iTemp)>ncol(out)) {
        out <- cbind(out,rep(NA,length(x)))
      }
      out[i,1:length(iTemp)] <- iTemp
    }
    names(out) <- paste0('V',c(1:ncol(out)))
    return(out)
  }
  
  #Function to parse out station ID and timeseries from prt files
  prtParse <- function(x) {
    idLine <- x[grep("Station -",x)[1]]
    Station.ID <- sub(".*?Station - (.*?) .*", "\\1", idLine)
    TS <- grep(pattern='^\\s{4}.[0-9]{4}\\s',x,value=T)
    return(list(Station.ID=Station.ID,
                TS = TS)
    )
  }
  
  
  #Parse out timeseries
  ##Get the text for each site
  textLines <- lapply(allFilesPRT,readLines)
  
  ##Parse out teh timeseries for each file
  TSData <- lapply(textLines,prtParse)
  
  ##Reformats the TS character vector into a dataframe
  siteTS <- lapply(TSData,tsFormat)
  
  ##Get site IDs for each file
  PRT_SiteID <- unlist(
    lapply(TSData,
           function(x) {
             return(x[1])
           })
  )
  ##rename lists with site IDs
  names(siteTS) <- paste0("X",PRT_SiteID)
  
  
  #subset GIS file to only sites found in EXP files
  
  ###Send warning message for dropped sites
  ###Get sites that are not in GIS file for warning message
  droppedSites <- gisData$Station.ID[which(!(gisData$Station.ID %in% PRT_SiteID))]
  
  if(length(droppedSites > 0))
  {
    warning(paste("The following sites are present in the GIS file and not found in PRT files in the Peak FQ output directory:",
                  droppedSites,sep=",")
    )
  }
  
  ###Subset GIS file to only sites that have an EXP and PRT
  gisData <- gisData[gisData$Station.ID %in% EXP_SiteID[,2],]
  gisData <- gisData[gisData$Station.ID %in% PRT_SiteID,]
  
  #This needs a comment
  recLen <- recCor <- matrix(NA,ncol=length(siteTS),nrow=length(siteTS))
  
  for (i in 1:length(siteTS)) {
    recLen[i,i] <- nrow(siteTS[[i]])
    idata <- abs(siteTS[[i]][,1:2])
    for (j in 1:i) {
      jdata <- abs(siteTS[[j]][,1:2])
      overlap <- intersect(idata[,1],jdata[,1])
      recLen[i,j] <- recLen[j,i] <- length(overlap)
      if (length(overlap)==0) {next}
      ijdata <- idata[which(is.element(idata[,1],overlap)),2]
      jidata <- jdata[which(is.element(jdata[,1],overlap)),2]
      if (length(unique(ijdata))==1|length(unique(jidata))==1) {next}
      recCor[i,j] <- recCor[j,i] <- cor(ijdata,jidata)
    }
  }
  
  #Make dataframes nice with column names, station IDs, and removing extraneous variables
  X <- gisData[,!names(gisData) %in% c("Lat","Long")]
  
  EXP_SiteID <- EXP_SiteID[2]
  colnames(EXP_SiteID) <- "Station.ID"
  
  AEP <- AEP[1,]
  
  colnames(Y) <- paste("AEP",AEP[1,],sep="_")
  Y$Station.ID <- EXP_SiteID$Station.ID
  Y <- Y[c(ncol(Y),1:ncol(Y)-1)]
  
  LP3f$Station.ID <- EXP_SiteID$Station.ID
  LP3f <- LP3f[c(ncol(LP3f),1:ncol(LP3f)-1)]
  
  colnames(LP3k) <- paste("AEP",AEP[1,],sep="_")
  LP3k$Station.ID <- EXP_SiteID$Station.ID
  LP3k <- LP3k[c(ncol(LP3k),1:ncol(LP3k)-1)]
  
  row.names(recLen) <- EXP_SiteID$Station.ID
  colnames(recLen) <- EXP_SiteID$Station.ID
  
  row.names(recCor) <- EXP_SiteID$Station.ID
  colnames(recCor) <- EXP_SiteID$Station.ID
  
  result <- list(
    sites=EXP_SiteID$Station.ID,
    Y=Y,
    AEP=AEP,
    X=X,
    LP3f=LP3f,
    LP3k=LP3k,
    BasChars=BasChars,
    MSEGR=unique(EXP_MSEGR[,2]),
    recLen=recLen,
    recCor=recCor
  )
  return(result)
  
}
