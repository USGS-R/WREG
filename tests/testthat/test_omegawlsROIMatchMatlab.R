context("Omega.WLS.ROImatchMatLab tests")

test_that("Run Omega.WLS.ROImatchMatLab",
          {
            # Import some example data
            expect_silent(
              
              load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
            )
            
            expect_silent(
              
              load(paste0(system.file("testData", package = "WREG"),"/omega.wls.roimatchmatlab.out.rda"))
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
            
            # Compute weighting matrix
            weightingResult <- Omega.WLS.ROImatchMatLab(Y.all = Y, X.all = X,
                                                        LP3.all = lp3Data, RecordLengths.all = staticData_peakFQ$recLen, NDX = NDX)
            
            expect_equal(weightingResult,omega.wls.roimatchmatlab.out)
            
            
            
          })



# test_that("Run Omega.GLS.ROImatchMatLab error handling",
#           {
#             # Import some example data
#             expect_silent(
#               
#               load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
#             )
#             
#             # Organizing input data
#             lp3Data <- staticData_peakFQ$LP3f
#             lp3Data$K <- staticData_peakFQ$LP3k$AEP_0.5
#             Y <- staticData_peakFQ$Y$AEP_0.5
#             X <- staticData_peakFQ$X[c("Sand", "OutletElev", "Slope")]
#             
#             #### Geographic Region-of-Influence
#             i <- 1 # Site of interest
#             n <- 10 # size of region of influence
#             Gdist <- vector(length=length(Y)) # Empty vector for geographic distances
#             for (j in 1:length(Y)) {
#               if (i!=j) {
#                 #### Geographic distance
#                 Gdist[j] <- Dist.WREG(Lat1 = staticData_peakFQ$BasChars$Lat[i],
#                                       Long1 = staticData_peakFQ$BasChars$Long[i],
#                                       Lat2 = staticData_peakFQ$BasChars$Lat[j],
#                                       Long2 = staticData_peakFQ$BasChars$Long[j]) # Intersite distance, miles
#               } else {
#                 Gdist[j] <- Inf # To block self identification.
#               }
#             }
#             temp <- sort.int(Gdist,index.return=TRUE)
#             NDX <- temp$ix[1:n] # Sites to use in this regression
#             
#             # Compute weighting matrix
#             expect_error(Omega.WLS.ROImatchMatLab(Y.all = "jazandapus", X.all = X,
#                                                         LP3.all = lp3Data, RecordLengths.all = staticData_peakFQ$recLen, NDX = NDX),
#                          "Invalid inputs were provided.  See warnings()."
#             )
#             expect_error(Omega.WLS.ROImatchMatLab(Y.all = Y , X.all = "jazandapus",
#                                                   LP3.all = lp3Data, RecordLengths.all = staticData_peakFQ$recLen, NDX = NDX),
#                          "missing value where TRUE/FALSE needed"
#             )
#             expect_error(Omega.WLS.ROImatchMatLab(Y.all = Y, X.all = X,
#                                                   LP3.all = "jazandapus", RecordLengths.all = staticData_peakFQ$recLen, NDX = NDX),
#                          "Invalid inputs were provided.  See warnings()."
#             )
#             expect_error(Omega.WLS.ROImatchMatLab(Y.all = Y, X.all = X,
#                                                   LP3.all = lp3Data, RecordLengths.all = "jazandapus", NDX = NDX),
#                          "argument is of length zero"
#             )
#             expect_error(Omega.WLS.ROImatchMatLab(Y.all = Y, X.all = X,
#                                                   LP3.all = lp3Data, RecordLengths.all = staticData_peakFQ$recLen, NDX = "jazandapus"),
#                          "Invalid inputs were provided.  See warnings()."
#             )
# 
#             
#             
#           })