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
#' \item{Y}{A matrix whose comlumns represent unique frequency events, 
#' while the row represent particular sites in the same order as \code{sites}.}
#' \item{X}{A matrix whose columns represent basin characteristics to be used
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
  siteInfo <- read.table(siteInfoFile,sep='\t',header=T)
  BasChars <- siteInfo[,is.element(names(siteInfo),
    c('Station.ID','Lat','Long'))]
  TestRecLen <- siteInfo$No..Annual.Series
  X <- siteInfo[,9:ncol(siteInfo)]
  
  # Load and parse FlowChar.txt
  flowCharFile <- file.path(wregPath,'FlowChar.txt')
  if (!file.exists(flowCharFile)) {
    stop(paste('Could not find',flowCharFile))
  }
  flowChar <- read.table(flowCharFile,sep='\t',header=T)
  if (!identical(siteInfo$Station.ID,flowChar$Station.ID)) {
    stop(paste0('The site order in ',siteInfoFile,' and ',
      flowCharFile,' must be identical.  ',
      '(All files must have identical order.)'))
  }
  Y <- flowChar[,2:ncol(flowChar)]
  names(Y) <- unlist(strsplit(names(Y),split='[.]'))
  
  # Load and parse LP3*.txt
  lp3gFile <- file.path(wregPath,'LP3G.txt')
  if (!file.exists(lp3gFile)) {
    stop(paste('Could not find',lp3gFile))
  }
  lp3g <- read.table(lp3gFile,sep='\t',header=T)
  if (!identical(siteInfo$Station.ID,lp3g$Station.ID)) {
    stop(paste0('The site order in ',siteInfoFile,' and ',
      lp3gFile,' must be identical.  ',
      '(All files must have identical order.)'))
  }
  test <- sum(unlist(
    lapply(lapply(t(lp3g[,2:ncol(lp3g)]),FUN=unique),FUN=length))>1)
  if (test>0) {
    stop(paste0('For a given set of tiem series, ',
      'the columns in LP3G.txt must be identical.'))
  }
  lp3kFile <- file.path(wregPath,'LP3K.txt')
  if (!file.exists(lp3kFile)) {
    stop(paste('Could not find',lp3kFile))
  }
  lp3k <- read.table(lp3kFile,sep='\t',header=T)
  if (!identical(siteInfo$Station.ID,lp3k$Station.ID)) {
    stop(paste0('The site order in ',siteInfoFile,' and ',
      lp3kFile,' must be identical.  ',
      '(All files must have identical order.)'))
  }
  lp3sFile <- file.path(wregPath,'LP3s.txt')
  if (!file.exists(lp3sFile)) {
    stop(paste('Could not find',lp3sFile))
  }
  lp3s <- read.table(lp3sFile,sep='\t',header=T)
  if (!identical(siteInfo$Station.ID,lp3s$Station.ID)) {
    stop(paste0('The site order in ',siteInfoFile,' and ',
      lp3sFile,' must be identical.  ',
      '(All files must have identical order.)'))
  }
  test <- sum(unlist(
    lapply(lapply(t(lp3s[,2:ncol(lp3s)]),FUN=unique),FUN=length))>1)
  if (test>0) {
    stop(paste0('For a given set of tiem series, ',
      'the columns in LP3s.txt must be identical.'))
  }
  LP3f <- data.frame(S=lp3s[,2],G=lp3g[,2],GR=siteInfo$Regional.Skew)
  LP3k <- lp3k[,2:ncol(lp3k)]
  
  # Screen for particular sites
  n <- nrow(siteInfo)
  ndx <- 1:nrow(siteInfo)
  if (length(sites)<2) {
    sites <- siteInfo$Station.ID
  }
  ndx <- which(is.element(siteInfo$Station.ID,as.integer(sites)))
  n2 <- length(ndx)
  if (n2!=length(sites)) {
    stop('You have requested sites that are not in SiteInfo.txt.')
  }
  sitesOut <- siteInfo$Station.ID[ndx]
  Y <- Y[ndx,]
  X <- X[ndx,]
  LP3f <- LP3f[ndx,]
  LP3k <- LP3k[ndx,]
  BasChars <- BasChars[ndx,]
  TestRecLen <- TestRecLen[ndx]
  
  # Load and parse UserWLS.txt (if available)
  uwlsFile <- file.path(wregPath,'UserWLS.txt')
  uwls <- NULL
  if (file.exists(uwlsFile)) {
    uwls <- read.table(uwlsFile,sep='\t',header=F)
    if (nrow(uwls)!=ncol(uwls)) {
      stop('UserWLS.txt must be a square matrix.')
    }
    if (nrow(uwls)!=n2) {
      uwls <- uwls[ndx,ndx]
      warning(paste0('UserWLS.txt was ammended to match sites in arguments.',
        ' (Note: Be sure weights do not depend on sites that were removed.)'))
    }
    uwls <- as.matrix(uwls)
    row.names(uwls) <- colnames(uwls) <- NULL
  }
  
  # Load and parse USGS*.txt files
  siteTS <- list()
  allFiles <- unique(c(apply(as.matrix(paste0('USGS0',sites,'*.txt')),
    MARGIN=1,FUN=list.files,path=wregPath),
    apply(as.matrix(paste0('USGS',sites,'*.txt')),
      MARGIN=1,FUN=list.files,path=wregPath)))
  if (length(allFiles)!=n2) {
    stop(paste0('Not all sites are represented with time series or ',
      'there are duplicate time series.'))
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
  if (!identical(diag(recLen),TestRecLen)) {
    stop(paste('Reported record lengths in SiteInfo.txt do not match',
      'time series files.'))
  }
  
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
