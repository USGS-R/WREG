---
title: "Weighted-Multiple-Linear Regression (WREG, v. 1.05-02)"
author: "William H. Farmer"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Weighted-Multiple-Linear Regression Program (WREG, v. 1.05-02)}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## 1 Introduction
The Weighted-Multiple-Linear Regression (WREG) package has been developed to replicate the functionality of the U.S. Geological Survey's standalone WREG program (1).  The WREG program was deployed as a MatLab executable program and is operated exclusively through a graphical user interface.  The WREG package is intended to provide an open-source resource for users wishing to apply the methods of the WREG program in a more accessible, adaptable and transparent environment using the R programming language.  The package is provided for users with some familiarity with the R programming environment.  For users who prefer a graphical user interface, this will be provided in a future release.

The WREG package and program are designed to help users develop regional regression equations relating streamflow frequency statistics to basin characteristics.  In addition to ordinary least-squares, WREG implements hydrologic applications of weighted and generalized least-squares (2 and 3).  The original manual for the WREG program (1) provides the necessary background and documentation of the WREG program and package.  This vignette guides potential users through the basic workflows required to implement ordinary, weighted and genearlized least-squares regression.

### 1.1 Compatability between the WREG program and package
In the development of the WREG package, several improvements on the WREG program were considered and implemented.  Full documentation of improvements will be included in a revised manual.  However, an option, `Legacy=TRUE` has been included to force the WREG package to return the same result as the WREG program.  This is only include for testing, but may help transitioning users become comfortable with the validity of the new implementations.

## 2 Package and example data
Once installed, the WREG package is loaded in the standard way:
```{r}
library(WREG)
```
The package contains the same sample input files that are included with the deployment of version 1.05 of the WREG program.  The input files have been condensed into list objects for easy manipulation and storage.  The data can be loaded with the following command:
```{r}
data(baseData)
```
The list object `baseData` contains the streamflow frequency statistics (oringally contained in FlowChar.txt), the independent variables (SiteInfo.txt) and the parameters of the fitted Log-Pearson Type III distributions (LP3G.txt, LP3K.txt, and LP3s.txt).

This example data set contains several dependent variables
```{r echo=FALSE}
names(baseData$Dependent[2:4])
```
and several independent variables
```{r echo=FALSE}
names(baseData$Independent[9:14])
```
allowing users to explore additional regressions and transformations beyond the scope of this vignette.

The `baseData` list also contains a matrix of the record lengths at each site and concurrent record lengths.  The WREG program derives these from analyzing the USGS########.txt files, but those have been pre-processed here.  (The orginal site files can be accessed with `data(siteData)`.)

## 3 Ordinary least-squares
Being the simplest form of regression implemented in the WREG package, ordinary least-squares regression reuires the fewest inputs into the main WREG function `WREG.MLR`.  Only independet and dependent variables are required to implement ordinary least-squares regression.

For this vignette we will consider a regression of the 1% peak streamflow against drainage area and precipitation.  We will fit a power-law regression, taking the natural logarithms of all independent and dependent variables.  The following line of code prepares the dependent variable.
```{r}
exY <- log(baseData$Dependent$Q1.)
```
The independent variables are defined as
```{r}
X1 <- log(baseData$Independent$DRNAREA)
X2 <- log(baseData$Independent$PRECIP)
X0 <- rep(1,length(X1))
exX <- cbind(X0,X1,X2)
```
Note that it is necessary to include the unit predictor in order to fit the leading coefficient or constant.

The WREG program provided a graphical user interface for the selection and transformation of variables.  This task is executed at the command line before passing inputs to the WREG package's main computational engine `WREG.MLR`.  Once the dependent and independent varaibles are defined, they are merely passed to `WREG.MLR` with the `Reg` option set to `'OLS'`:
```{r}
Ex.OLS <- WREG.MLR(Y=exY,X=exX,Reg='OLS')
```
The resulting list object, `Ex.OLS` contains a wide array of outputs
```{r echo=FALSE}
names(Ex.OLS)
```
each of which is described in detail in the help file of `WREG.MLR`.  The outputs of the WREG package are designed to be more comprehensive than the outputs of the WREG program.

The coefficients, and associated statistics thereof, are included in `Coefs`
```{r echo=FALSE}
knitr::kable(round(Ex.OLS$Coefs,digits=4))
```

Performance metrics are found in `PerformanceMetrics`
```{r echo=FALSE}
knitr::kable(round(data.frame(t(unlist(Ex.OLS$PerformanceMetrics))),digits=4))
```

Though not presented here, the residuals, leverage and inlufence are found in `ResLevInf`.

These outputs can be used to create a wide range of performance plots for detailed analysis.  Here are several plots mirroring the output of the WREG program.  (It is important to remember that the units of the axes replicate the units of the dependent variable `Y`.)

```{r eval=FALSE, fig.show='hold', fig.width = 6.5, fig.height = 6.5,fig.align='center'}
plot(Ex.OLS$fitted.values,Ex.OLS$residuals,
  main='Residuals versus Estimated Flow Characteristics',
  xlab='Estimated Flow Characteristic',ylab='Residual',
  xaxs='i',yaxs='i',ylim=c(-1.5,1.5),xlim=c(4.5,10))
grid(nx=11,ny=6)
```

```{r eval=FALSE, fig.show='hold', fig.width = 6.5, fig.height = 6.5,fig.align='center'}
plot(Ex.OLS$Y,Ex.OLS$ResLevInf$Leverage,
  main='Leverage Values versus Observations',
  xlab='Observation',ylab='Leverage Value',
  xaxs='i',yaxs='i',ylim=c(0,0.3),xlim=c(5,11))
grid(nx=6,ny=6)
abline(Ex.OLS$LevLim,0,col='red')
```

```{r eval=FALSE, fig.show='hold', fig.width = 6.5, fig.height = 6.5,fig.align='center'}
plot(Ex.OLS$Y,Ex.OLS$ResLevInf$Influence,
  main='Influence Value versus Observation',
  xlab='Observation',ylab='Influence Value',
  xaxs='i',yaxs='i',ylim=c(0,0.14),xlim=c(5,11))
grid(nx=6,ny=7)
abline(Ex.OLS$InflLim,0,col='red')
```

## 4 Weighted Least-Squares
With a few additional inputs, `WREG.MLR` can apply weighted least-squares regression as defined by (2).  In order to compute weights, the function requires information about record lengths and the fitted Log-Pearson Type III distribution.

Example record lengths are included in the example data set:
```{r}
RL <- baseData$RecordLengths
```

The example record lengths include the concurrent record lengths, resulting in a matrix.  A vector of record lengths, ignoring the concurrent records, is obtained with `diag(RL)`; if a matrix is passed to `WREG.MLR` when `Reg='WLS'`, requesting weighted least-squares, then the diagnoal is extracted automatically.

In addition to record lengths, the weighting designed by (2) requires information on the fitted Log-Pearson Type III distribution.  The three parameters, the standard deviation, deviates and skew, are passed as a list to `WREG.MLR`.  The example files contain the LP3 parameters for each dependent variable, so the following code extracts those realted to the 1% event.
```{r}
exLP3 <- list(S=baseData$LP3$S$s.1.,K=baseData$LP3$K$K.50.,G=baseData$LP3$G$Skew50.)
```

With these inputs, it is then possible to compute the coefficients based on the weight least-squares algorithm:
```{r}
Ex.WLS <- WREG.MLR(Y=exY,X=exX,Reg='WLS',RecordLengths=RL,LP3=exLP3)
```

The new coefficients are
```{r echo=FALSE}
knitr::kable(round(Ex.WLS$Coefs,digits=4))
```

The application of weighted least-squared regression produces a range of new performance metrics,
```{r echo=FALSE}
names(Ex.WLS$PerformanceMetrics)
```

## 5 Generalized Least-Squares
The final form of regression implemented in the WREG package is the hydrologic interpretation of generalized least-squares regression developed by (3).  This application of generalized least-squares regression relies on an approximation of cross correlation based on inter-site distances.  With this approximation, generalized least-squares regression can be applied with or without uncertainty in regional skew.

First we need to collect station IDs and latitude-longitude couplets.
```{r}
descBasin <- baseData$Independent[,1:3]
names(descBasin)
```

These basin characteristics are used to either assess the fit of estimated cross-correlation or apply generalized least-squares regression.

### 5.1 Estimated cross-correlation
The WREG program provides a graphical user interface in order to determine a relationship between distance and correlation.  In the WREG package, the function `testXCorr` can be used to undertake the same exploration.  Site data, like that provided in `siteData` are needed to assess model fits.

```{r}
data(siteData)
```

By manually varying the `alpha`, `theta` and `concurrentMin`, one can assess the fitted estimates.  The model can also be fit using two separate calculations of distances: the WREG program uses a nautical-mile approximation, while 'DistMeth=2' applies the Haversine formula.

```{r fig.show='hold', fig.width = 6.5, fig.height = 6.5,fig.align='center'}
testXCorr(siteData=siteData,BasinChars=descBasin,alpha=0.002,theta=0.98,DistMeth=2,
  concurrentMin=25,plot=TRUE)
```

### 5.2 Without uncertainty in skew
Once `alpha` and `theta` have been selected, it is possible to fit a model using genearlized least-squares.

```{r}
Ex.GLS <- WREG.MLR(Y=exY,X=exX,Reg='GLS',RecordLengths=RL,LP3=exLP3,
                     BasinChars=descBasin,alpha=0.002,theta=0.98,DistMeth=2)
```

The new coefficients are
```{r echo=FALSE}
knitr::kable(round(Ex.GLS$Coefs,digits=4))
```

The same perfromance metrics calculated when appolying wegihted least-squares regression are calculated when computing generalized least-squares regression.

### 5.3 Uncertainty in skew
In addition to the generalized least-squares developed by (3), the WREG package is capable of implementing the correction for regional uncertainty developed by (4).  When requesting the skew-correction, `Reg='GLSskew'`, it is necessary to include several additional parameters:

- `MSEGR` : The mean-squared-error of the regional skew
- `TY` : The return period of the streamflow frequency statistic
- `Peak` : A logical value indicating if the event is a peak event (`TRUE`) or a low event (`FALSE`)

It is also necessary to add the regional skew value to the LP3 list:

```{r}
exLP3$GR <- baseData$Independent$Regional.Skew
```

Applying generalized least-squares regression with a correction for uncertainty in regional skew

```{r}
Ex.GLSs <- WREG.MLR(Y=exY,X=exX,Reg='GLSskew',RecordLengths=RL,LP3=exLP3,
  BasinChars=descBasin,alpha=0.002,theta=0.98,DistMeth=2,
  MSEGR=0.302,TY=100,Peak=TRUE)
```

produces a new set of coefficients

```{r echo=FALSE}
knitr::kable(round(Ex.GLSs$Coefs,digits=4))
```

## 6 Region-of-influence regression
The WREG program is also capable of applying region-of-influence regression, calculating coefficients by any of the methods previously discussed.  This same functionality is replicated in the WREG package using the `WREG.RoI` function.  This vignette will not explore the application of region-of-influence regression, but the associated help files may provide the guidance needed.  Naturally, applying region-of-influence produces several unique outputs; again, the help files provide descriptions.

## 7 References

1. Eng, K., Chen, Y., and Kiang, J.E., 2009, User's guide to the weighted-multiple-linear-regression program (WREG version 1.0): U.S. Geological Survey Techniques and Methods, book 4, chap. A8, 21 p.
1. Tasker, G.D., 1980, Hydrologic regression with weighted least squares: Water Resources Research, v. 16, no. 6, p. 1107-1113.
1. Stedinger, J.R, and Tasker, G.D., 1985, Regional hydrologic analysis, 1, ordinary, weighted, and generalized least squares compared: Water Resources Research, v. 21, no. 9, p. 1421-1432.
1. Griffis, V.W., and Stedinger, J.R., 2007, The use of GLS regression in regional hydrologic analyses: Journal of Hydrology, v. 344, p. 82-95.