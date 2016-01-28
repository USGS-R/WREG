## ------------------------------------------------------------------------
library(WREG)

## ------------------------------------------------------------------------
data(baseData)

## ----echo=FALSE----------------------------------------------------------
names(baseData$Dependent[2:4])

## ----echo=FALSE----------------------------------------------------------
names(baseData$Independent[9:14])

## ------------------------------------------------------------------------
exY <- log(baseData$Dependent$Q1.)

## ------------------------------------------------------------------------
X1 <- log(baseData$Independent$DRNAREA)
X2 <- log(baseData$Independent$PRECIP)
X0 <- rep(1,length(X1))
exX <- cbind(X0,X1,X2)

## ------------------------------------------------------------------------
Ex.OLS <- WREG.MLR(Y=exY,X=exX,Reg='OLS')

## ----echo=FALSE----------------------------------------------------------
names(Ex.OLS)

## ----echo=FALSE----------------------------------------------------------
knitr::kable(round(Ex.OLS$Coefs,digits=4))

## ----echo=FALSE----------------------------------------------------------
knitr::kable(round(data.frame(t(unlist(Ex.OLS$PerformanceMetrics))),digits=4))

## ----eval=FALSE, fig.show='hold', fig.width = 6.5, fig.height = 6.5,fig.align='center'----
#  plot(Ex.OLS$fitted.values,Ex.OLS$residuals,
#    main='Residuals versus Estimated Flow Characteristics',
#    xlab='Estimated Flow Characteristic',ylab='Residual',
#    xaxs='i',yaxs='i',ylim=c(-1.5,1.5),xlim=c(4.5,10))
#  grid(nx=11,ny=6)

## ----eval=FALSE, fig.show='hold', fig.width = 6.5, fig.height = 6.5,fig.align='center'----
#  plot(Ex.OLS$Y,Ex.OLS$ResLevInf$Leverage,
#    main='Leverage Values versus Observations',
#    xlab='Observation',ylab='Leverage Value',
#    xaxs='i',yaxs='i',ylim=c(0,0.3),xlim=c(5,11))
#  grid(nx=6,ny=6)
#  abline(Ex.OLS$LevLim,0,col='red')

## ----eval=FALSE, fig.show='hold', fig.width = 6.5, fig.height = 6.5,fig.align='center'----
#  plot(Ex.OLS$Y,Ex.OLS$ResLevInf$Influence,
#    main='Influence Value versus Observation',
#    xlab='Observation',ylab='Influence Value',
#    xaxs='i',yaxs='i',ylim=c(0,0.14),xlim=c(5,11))
#  grid(nx=6,ny=7)
#  abline(Ex.OLS$InflLim,0,col='red')

## ------------------------------------------------------------------------
RL <- baseData$RecordLengths

## ------------------------------------------------------------------------
exLP3 <- list(S=baseData$LP3$S$s.1.,K=baseData$LP3$K$K.50.,G=baseData$LP3$G$Skew50.)

## ------------------------------------------------------------------------
Ex.WLS <- WREG.MLR(Y=exY,X=exX,Reg='WLS',RecordLengths=RL,LP3=exLP3)

## ----echo=FALSE----------------------------------------------------------
knitr::kable(round(Ex.WLS$Coefs,digits=4))

## ----echo=FALSE----------------------------------------------------------
names(Ex.WLS$PerformanceMetrics)

## ------------------------------------------------------------------------
descBasin <- baseData$Independent[,1:3]
names(descBasin)

## ------------------------------------------------------------------------
data(siteData)

## ----fig.show='hold', fig.width = 6.5, fig.height = 6.5,fig.align='center'----
testXCorr(siteData=siteData,BasinChars=descBasin,alpha=0.002,theta=0.98,DistMeth=2,
  concurrentMin=25,plot=TRUE)

## ------------------------------------------------------------------------
Ex.GLS <- WREG.MLR(Y=exY,X=exX,Reg='GLS',RecordLengths=RL,LP3=exLP3,
                     BasinChars=descBasin,alpha=0.002,theta=0.98,DistMeth=2)

## ----echo=FALSE----------------------------------------------------------
knitr::kable(round(Ex.GLS$Coefs,digits=4))

## ------------------------------------------------------------------------
exLP3$GR <- baseData$Independent$Regional.Skew

## ------------------------------------------------------------------------
Ex.GLSs <- WREG.MLR(Y=exY,X=exX,Reg='GLSskew',RecordLengths=RL,LP3=exLP3,
  BasinChars=descBasin,alpha=0.002,theta=0.98,DistMeth=2,
  MSEGR=0.302,TY=100,Peak=TRUE)

## ----echo=FALSE----------------------------------------------------------
knitr::kable(round(Ex.GLSs$Coefs,digits=4))

