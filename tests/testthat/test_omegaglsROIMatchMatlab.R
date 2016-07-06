context("Omega.GLS.ROImatchMatLab tests")

test_that("Run Omega.GLS.ROImatchMatLab",
          {
            # Import some example data
            expect_silent(
              
              load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
            )
            
            expect_silent(
              
              load(paste0(system.file("testData", package = "WREG"),"/omega.gls.roimatchmatlab.out.rda"))
            )
            # Organizing input data
            lp3Data <- staticData_peakFQ$LP3f
            lp3Data$K <- staticData_peakFQ$LP3k$AEP_0.5
            Y <- staticData_peakFQ$Y$AEP_0.5
            X <- staticData_peakFQ$X[c("Sand", "OutletElev", "Slope")]
            
            #### Geographic Region-of-Influence
            i <- 1 # Site of interest
            n <- 10 # size of region of influence
            Gdist <- vector(length=length(Y)) # Empty vector for geographic distances
            for (j in 1:length(Y)) {
              if (i!=j) {
                #### Geographic distance
                Gdist[j] <- Dist.WREG(Lat1 = staticData_peakFQ$BasChars$Lat[i],
                                      Long1 = staticData_peakFQ$BasChars$Long[i],
                                      Lat2 = staticData_peakFQ$BasChars$Lat[j],
                                      Long2 = staticData_peakFQ$BasChars$Long[j]) # Intersite distance, miles
              } else {
                Gdist[j] <- Inf # To block self identification.
              }
            }
            temp <- sort.int(Gdist,index.return=TRUE)
            NDX <- temp$ix[1:n] # Sites to use in this regression
            
            # Pull out characeristics of the region.
            Y.i <- Y[NDX] # Predictands from region of influence
            X.i <- X[NDX,] # Predictors from region of influence
            RecordLengths.i <- staticData_peakFQ$recLen[NDX,NDX] # Record lengths from region of influence
            BasinChars.i <- staticData_peakFQ$BasChars[NDX,] # Basin characteristics (IDs, Lat, Long) from region of influence.
            LP3.i <- data.frame(lp3Data)[NDX,] # LP3 parameters from region of influence
            
            # Compute weighting matrix
            weightingResult <- Omega.GLS.ROImatchMatLab(alpha = 0.01, theta = 0.98,
                                                        Independent = staticData_peakFQ$BasChars, X = X.i, Y = Y.i,
                                                        RecordLengths = RecordLengths.i, LP3 = LP3.i, TY = 20,
                                                        X.all = X, LP3.all = lp3Data)
            expect_equal(weightingResult,omega.gls.roimatchmatlab.out)
            
            
            
          })
