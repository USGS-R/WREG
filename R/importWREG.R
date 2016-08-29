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
importWREG <- function(wregPath,sites='') {
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
  siteInfo <- read.table(siteInfoFile,sep='\t',header=T,
                         colClasses = list(Station.ID='character'))
  siteInfo$Station.ID <- ifelse(nchar(siteInfo$Station.ID)%%2>0,
                                paste0("0",siteInfo$Station.ID),siteInfo$Station.ID)
  BasChars <- siteInfo[,is.element(names(siteInfo),
                                   c('Station.ID','Lat','Long'))]
  X <- siteInfo[,c(1,9:ncol(siteInfo))]
  
  # Screen for particular sites
  n <- nrow(siteInfo)
  ndx <- 1:nrow(siteInfo)
  if (length(sites)<2) {
    sites <- siteInfo$Station.ID
  }
  ndx <- which(is.element(siteInfo$Station.ID,sites))
  n2 <- length(ndx)
  if (n2!=length(sites)) {
    stop('You have requested sites that are not in SiteInfo.txt.')
  }
  sitesOut <- siteInfo$Station.ID[ndx]
  X <- X[ndx,]
  BasChars <- BasChars[ndx,]
  TestRecLen <- siteInfo[,is.element(names(siteInfo),
                                     c('Station.ID','No..Annual.Series'))]
  TestRecLen <- TestRecLen[ndx,]
  site1 <- sitesOut
  
  # Load and parse FlowChar.txt
  flowCharFile <- file.path(wregPath,'FlowChar.txt')
  if (!file.exists(flowCharFile)) {
    stop(paste('Could not find',flowCharFile))
  }
  Y <- read.table(flowCharFile,sep='\t',header=T,
                  colClasses = list(Station.ID='character'))
  # Y <- Y[sort.int(runif(nrow(Y)),index.return = TRUE)$ix,] # for testing
  Y$Station.ID <- ifelse(nchar(Y$Station.ID)%%2>0,
                         paste0("0",Y$Station.ID),Y$Station.ID)
  names(Y)[2:ncol(Y)] <- unlist(strsplit(names(Y)[2:ncol(Y)],split='[.]'))
  ndx <- which(is.element(Y$Station.ID,sitesOut))
  Y <- Y[ndx,]
  site2 <- Y$Station.ID
  Y <- Y[match(sitesOut,Y$Station.ID),]
  
  # Load and parse LP3*.txt
  lp3gFile <- file.path(wregPath,'LP3G.txt')
  if (!file.exists(lp3gFile)) {
    stop(paste('Could not find',lp3gFile))
  }
  lp3g <- read.table(lp3gFile,sep='\t',header=T)
  lp3g$Station.ID <- ifelse(nchar(lp3g$Station.ID)%%2>0,
                            paste0("0",lp3g$Station.ID),lp3g$Station.ID)
  test <- sum(unlist(
    lapply(lapply(t(lp3g[,2:ncol(lp3g)]),FUN=unique),FUN=length))>1)
  if (test>0) {
    stop(paste0('For a given set of tiem series, ',
                'the columns in LP3G.txt must be identical.'))
  }
  ndx <- which(is.element(lp3g$Station.ID,sitesOut))
  lp3g <- lp3g[ndx,]
  site3a <- lp3g$Station.ID
  lp3g <- lp3g[match(sitesOut,lp3g$Station.ID),]
  lp3kFile <- file.path(wregPath,'LP3K.txt')
  if (!file.exists(lp3kFile)) {
    stop(paste('Could not find',lp3kFile))
  }
  lp3k <- read.table(lp3kFile,sep='\t',header=T)
  lp3k$Station.ID <- ifelse(nchar(lp3k$Station.ID)%%2>0,
                            paste0("0",lp3k$Station.ID),lp3k$Station.ID)
  ndx <- which(is.element(lp3k$Station.ID,sitesOut))
  lp3k <- lp3k[ndx,]
  site3b <- lp3k$Station.ID
  lp3k <- lp3k[match(sitesOut,lp3k$Station.ID),]
  lp3sFile <- file.path(wregPath,'LP3S.txt')
  if (!file.exists(lp3sFile)) {
    stop(paste('Could not find',lp3sFile))
  }
  lp3s <- read.table(lp3sFile,sep='\t',header=T)
  lp3s$Station.ID <- ifelse(nchar(lp3s$Station.ID)%%2>0,
                            paste0("0",lp3s$Station.ID),lp3s$Station.ID)
  test <- sum(unlist(
    lapply(lapply(t(lp3s[,2:ncol(lp3s)]),FUN=unique),FUN=length))>1)
  if (test>0) {
    stop(paste0('For a given set of tiem series, ',
                'the columns in LP3S.txt must be identical.'))
  }
  ndx <- which(is.element(lp3s$Station.ID,sitesOut))
  lp3s <- lp3s[ndx,]
  site3c <- lp3s$Station.ID
  lp3s <- lp3s[match(sitesOut,lp3s$Station.ID),]
  
  LP3f <- data.frame(Station.ID=sitesOut,S=lp3s[,2],G=lp3g[,2],
                     GR=siteInfo$Regional.Skew[ndx])
  LP3k <- lp3k
  
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
  
  # Load and parse USGS*.txt files
  siteTS <- list()
  allFiles <- unique(apply(as.matrix(paste0('USGS',sitesOut,'*.txt')),
                           MARGIN=1,FUN=list.files,path=wregPath))
  if (any(lapply(allFiles,length) == 0)) {
    warning(paste("The following sites do not have timeseries data",
                  sitesOut[which(lapply(allFiles,length) == 0)]))
    #Remove missing file from all files list
    sitesOut <- sitesOut[-which(lapply(allFiles,length) == 0)]
    allFiles <- allFiles[-which(lapply(allFiles,length) == 0)]
    
    #Remove missing site from all files
    Y <- Y[Y$Station.ID %in% sitesOut,]
    X <- X[X$Station.ID %in% sitesOut,]
    LP3f <- LP3f[LP3f$Station.ID %in% sitesOut,]
    LP3k <- LP3k[LP3k$Station.ID %in% sitesOut,]
    BasChars <- BasChars[BasChars$Station.ID %in% sitesOut,]
    siteInfo <- siteInfo[siteInfo$Station.ID %in% sitesOut,]
    TestRecLen <- TestRecLen[TestRecLen$Station.ID %in% sitesOut,]
    
  }
  if(length(allFiles) > length(sitesOut))
  {
    stop("Duplicate timeseries were found.")
  }
  allFiles <- file.path(wregPath,allFiles)
  siteTS <- lapply(allFiles,read.table)
  
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
      recCor[i,j] <- recCor[j,i] <- cor(ijdata,jidata)
    }
  }
  tester <- diag(recLen)
  names(tester) <- NULL
  if (!identical(tester,TestRecLen[,2])) {
    stop(paste('Reported record lengths in SiteInfo.txt do not match',
               'time series files.'))
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
