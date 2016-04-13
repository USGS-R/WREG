
#'Exploring correlation function
#'
#'@description \code{xcorPlot} explores the fit of the proposed approximation of
#'  cross-correlation.
#'  
#'@param object A list of the imported data from \code{\link{importPeakFQ}} 
#'  or\code{\link{importWREG}}.
#'@param alpha A number, \code{alpha} is a parameter used in the estimated 
#'  cross-correlation between site records.  See equation 20 in the WREG v. 1.05
#'  manual.
#'@param theta A number, \code{theta} is a parameter used in the estimated 
#'  cross-correlation between site records.  See equation 20 in the WREG v. 1.05
#'  manual.
#'@param concurrentMin A number specifying the minimum number of years of 
#'  concurrent record required to estimate cross-correlation.
#'@param DistMeth A value of \code{1} indicates that the "Nautical Mile" 
#'  approximation should be used to calculate inter-site distances.  A value of 
#'  \code{2} designates the Haversine approximation.  See 
#'  \code{\link{Dist.WREG}}.
#'@param plot A logical spcifying if the plot should be created.
#'  
#'@return If \code{plot=FALSE}, the Nash-Sutcliffe model efficieny is returned.
#'  
#'@import stats graphics
#'  
#' @examples 
#' \dontrun Add example
#'@export
xcorPlot <- function(object,alpha,theta,concurrentMin,
  DistMeth=2,plot=TRUE) {
  # William Farmer, October 22, 2015
  # Revised by WHF, March 02, 2016
  
  if (concurrentMin < 10) {
    warning(paste0('It is not reccommended to use a concurrent record',
      ' length less than 10 years. The value has been increased to 10.'))
    concurrentMin <- 10
  }
  
  n <- nrow(object$BasChars)
  plotData <- matrix(NA,ncol=2,nrow=n^2)
  iter <- 0
  maxDist <- maxcor <- -Inf
  mincor <- Inf
  for (i in 1:n) {
    for (j in i:n) {
      if (i!=j) {
        ijDist <- Dist.WREG(object$BasChars$Lat[i],object$BasChars$Long[i],
          object$BasChars$Lat[j],object$BasChars$Long[j],method=DistMeth)
        if (object$recLen[i,j]>9) {
          maxDist <- max(maxDist,ijDist)
        }
        if (object$recLen[i,j] >=concurrentMin) {
          iter <- iter + 1
          plotData[iter,1] <- ijDist
          plotData[iter,2] <- object$recCor[i,j]
        }
      }
    }
  }
  
  ndx <- which(!is.na(rowSums(plotData)))
  plotData <- plotData[ndx,]
  tester <- object$recCor
  tester[object$recLen<10] <- NA
  diag(tester) <- NA
  mincor <- floor(min(c(tester),na.rm=T)*10)/10
  maxcor <- ceiling(max(c(tester),na.rm=T)*10)/10
  
  estRhos <- theta^(plotData[,1]/(alpha*plotData[,1]+1))
  
  nse <- 1 - sum((estRhos-plotData[,2])^2)/
    sum((plotData[,2]-mean(plotData[,2]))^2)
  
  if (plot) {
    ny <- round((maxcor-mincor)/0.1)
    maxpower <- 10
    splits <- 10^0
    for (i in seq(0.2,1,0.1)) {
      splits <- c(splits,10^seq(log10(i*10^1),log10(i*10^maxpower),1))
    }
    splits <- sort(splits)
    splits1 <- splits-(maxDist)
    splits1[splits1<0] <- NA
    ndx <- which(splits1==min(splits1,na.rm=T))
    xlim2 <- round(splits[ndx])
    plotDists <- seq(0,xlim2,length.out=1000)
    plotRhos <- theta^(plotDists/(alpha*plotDists+1))
    xlim <- c(0,xlim2)
    graphics::plot(plotData[,1],plotData[,2],type='p',pch=1,
      main=paste0('Correlation Smoothing Function\n',
        '(alpha=',alpha,', theta=',theta,', NSE=',round(nse,digits=4),')'),
      xlab='Geographic distance (km)',
      ylab='Sample correlation',
      xaxs='i',yaxs='i',
      ylim=c(mincor,maxcor),xlim=xlim)
    graphics::lines(plotDists,plotRhos,lty=5,col='red')
    graphics::grid(nx=round(xlim2/10^floor(log10(xlim2))*4),ny=ny)
    graphics::legend('topright',
      legend=c(paste0('Observation (Y=',concurrentMin,')'),'Model'),lty=c(NA,5),
      col=c('black','red'),pch=c(1,NA))
  } else {
    return(nse)
  }
  
  
  
}