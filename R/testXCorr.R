
#' Exploring correlation function
#'
#'@description
#'\code{testXCorr} explores the fit of the proposed approximation of 
#'cross-correlation.
#'
#' @param siteData A list of the records at each site.  Each element of the list,
#' representing a site, is a data frame with 2 variables names \code{year} 
#' and \code{obs}.
#' @param BasinChars A dataframe containing three variables: \code{StationID} is the 
#' numerical identifier (without a leading zero) of each site, \code{Lat} is the latitude
#' of the site, in decimal degrees, and \code{Long} is the longitude of the site, in decimal
#' degrees.
#' @param alpha A number, \code{alpha} is a parameter used in the estimated 
#' cross-correlation between site records.  See equation 20 in the WREG 
#' v. 1.05 manual.
#' @param theta A number, \code{theta} is a parameter used in the estimated 
#' cross-correlation between site records.  See equation 20 in the WREG 
#' v. 1.05 manual.
#' @param DistMeth A value of \code{1} indicates that the "Nautical Mile" 
#' approximation should be used to calculate inter-site distances.  A value 
#' of \code{2} designates the Haversine approximation.  See 
#' \code{\link{Dist.WREG}}.
#' @param concurrentMin A number specifying the minimum number of years of 
#' concurrent record required to estimate cross-correlation.
#' @param plot A logical spcifying if the plot should be created.
#' 
#' @return If \code{plot=FALSE}, the Nash-Sutcliffe model efficieny is returned.
#' 
#' @import stats graphics
#' @export
testXCorr <- function(siteData,BasinChars,alpha,theta,DistMeth=c(1,2),
  concurrentMin=25,plot=TRUE) {
  # William Farmer, October 22, 2015
  #
  
  plotData <- matrix(NA,ncol=2,nrow=length(siteData)^2)
  RecordLengths <- matrix(NA,ncol=length(siteData),nrow=length(siteData))
  iter <- 0
  maxDist <- maxcor <- -Inf
  mincor <- Inf
  for (i in 1:length(siteData)) {
    iData <- siteData[[i]]
    for (j in i:length(siteData)) {
      if (i!=j) {
        ijDist <- Dist.WREG(BasinChars$Lat[i],BasinChars$Long[i],
          BasinChars$Lat[j],BasinChars$Long[j],method=DistMeth)
        maxDist <- max(maxDist,ijDist)
        jData <- siteData[[j]]
        ### Concurrent period of record
        interY <- intersect(iData$year,jData$year)
        RecordLengths[i,j] <- length(interY)
        if (RecordLengths[i,j]>2) {
          indx <- which(is.element(iData$year,interY))
          jndx <- which(is.element(jData$year,interY))
          ijcor <- NA
          if (sd(iData$obs[indx])>0&sd(jData$obs[jndx])>0) {
            ijcor <- stats::cor(iData$obs[indx],jData$obs[jndx],
              method='pearson')
          }
          mincor <- min(mincor,ijcor,na.rm=T)
          maxcor <- max(maxcor,ijcor,na.rm=T)
          if (RecordLengths[i,j]>=concurrentMin) {
            # if above minimum, add to plot data
            iter <- iter + 1
            plotData[iter,1] <- ijDist
            plotData[iter,2] <- ijcor
          }
        }
      }
    }
  }
  
  ndx <- which(!is.na(rowSums(plotData)))
  plotData <- plotData[ndx,]
  
  estRhos <- theta^(plotData[,1]/(alpha*plotData[,1]+1))
  
  nse <- 1 - sum((estRhos-plotData[,2])^2)/
    sum((plotData[,2]-mean(plotData[,2]))^2)
  
  if (plot) {
    ylim1 <- ifelse(round(mincor,digits=1)>mincor,round(mincor,digits=1)-0.1,
      round(mincor,digits=1))
    ylim2 <- ifelse(round(maxcor,digits=1)<maxcor,round(maxcor,digits=1)+0.1,
      round(maxcor,digits=1))
    ylim <- c(ylim1,ylim2)
    ny <- round((ylim2-ylim1)/0.1)
    nx <- 10
    splits <- c(1,5,10,50,100,500,1000,5000,10000)
    splits1 <- splits-(maxDist/nx)
    splits1[splits1<0] <- NA
    ndx <- which(splits1==min(splits1,na.rm=T))
    xlim2 <- splits[ndx]*nx
    plotDists <- seq(0,xlim2,length.out=1000)
    plotRhos <- theta^(plotDists/(alpha*plotDists+1))
    xlim <- c(0,xlim2)
    graphics::plot(plotData[,1],plotData[,2],type='p',
      main=paste0('Correlation Smoothing Function\n',
        '(alpha=',alpha,', theta=',theta,', NSE=',round(nse,digits=4),')'),
      xlab='Geographic Distance (km)',
      ylab='Sample rho',
      xaxs='i',yaxs='i',
      ylim=ylim,xlim=xlim)
    graphics::lines(plotDists,plotRhos,lty=5,col='red')
    graphics::grid(nx=nx,ny=ny)
    graphics::legend('topright',
      legend=c(paste0('Observation (Y=',concurrentMin,')'),'Model'),lty=c(NA,5),
      col=c('black','red'),pch=c(0,NA))
  } else {
    return(nse)
  }
  
  
  
}