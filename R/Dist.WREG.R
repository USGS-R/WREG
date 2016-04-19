#'Distance calculation (WREG)
#'
#'@description \code{Dist.WREG} calculates the distance between two points 
#'  defined by Latitude-Longitude coordinates in decimal degrees.
#'  
#'@param Lat1 Latitude of the first point, in decimal degrees.
#'@param Long1 Lonigtude of the first point, in decimal degrees.
#'@param Lat2 Latitude of the second point, in decimal degrees.
#'@param Long2 Longitude of the second point, in deceimal degrees.
#'@param method Idicates which technique to use for distance calculation.  See 
#'  details.
#'  
#'@details \code{Dist.WREG} is capable of using two techniques to calculate 
#'  intersite distances. \code{method==1} indicates that the "Nautical Mile" 
#'  approximation should be used. This is the function that is currently
#'  employed by WREG v. 1.05.  Each arcminute is equal to 1852 meters.
#'  \code{method==2} indicates that the Haversine approximation should be used.
#'  
#'@return Returns the distance between the two sites in miles.
#'  
#' @examples 
#' \dontrun{
#' #add examples
#' }
#'@export
Dist.WREG <- function(Lat1,Long1,Lat2,Long2,method=2) {
  # William Farmer, USGS, January 23, 2015
  
  # basic error checking
  ivars <- c('Lat1','Long1','Lat2','Long2')
  err <- FALSE
  for (i in ivars) {
    check <- eval(parse(text=paste0('missing(',i,')')))
    if (check) {
      err = TRUE
      warning(paste(i,"must be provided."))
      next
    }
    check <- eval(parse(text=paste0('sum(!is.numeric(',i,'))>0')))
    if (check) {
      err = TRUE
      warning(paste(i,"is not numeric."))
      next
    }
    check <- eval(parse(text=paste0('sum(is.infinite(',i,'))>0')))
    if (check) {
      err = TRUE
      warning(paste(i,"contains infinite values."))
    }
    check <- eval(parse(text=paste0('sum(is.na(',i,'))>0')))
    if (check) {
      err = TRUE
      warning(paste(i,"contains missing values."))
    }
  }
  if (!is.element(method,c(1,2))) {
    warning("'method' must be either 1 for use of a nautical mile ",
      "approximation or 2 for use of the haversine formula.")
    err <- TRUE
  }
  if (err) {
    stop("Invalid inputs were provided.  See warnings()")
  }
  
  if (method==1) {
    ## Nautical mile conversion
    Dist <- sqrt((abs(Lat2-Lat1)*1852*60)^2+(abs(Long2-Long1)*1852*60)^2)*0.6214/1000 # Intersite distance, miles
  } else if (method==2) {
    R <- 6371/1.609 # Radius of the earth in miles
    ## Convert Lat/Long to Radians
    Lat1 <- Lat1*pi/180
    Long1 <- Long1*pi/180
    Lat2 <- Lat2*pi/180
    Long2 <- Long2*pi/180
    ## Haversine formula
    a <- (sin(0.5*(Lat2-Lat1)))^2+cos(Lat1)*cos(Lat2)*(sin(0.5*(Long2-Long1))^2)
    c <- 2*atan2(sqrt(a),sqrt(1-a))
    Dist <- R*c # Intersite distance, miles
  }
  return(Dist)
}