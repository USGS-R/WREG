#' Utility function to add and remove warnings without halting execution
#' 
#' @description
#' The \code{warn} handles warnings
#' for WREG-R.
#' 
#' @param mode Whether to initialize, check, get, add, or clear warnings
#' @param w (optional) Message for adding warning
#' 
#' @details
#' This functions is a warning scheme that does not halt execution.
#' 
#' @return Nothing, hasWarnings, or string list of warnings
#' 
#'@examples
#'hasWarnings <- warn("check")
#'
#'@export

warn <- function(mode, w=NULL){

  init <- function(clear=FALSE){
    
    if(clear){
      wregWarnings <<- list()
      hasWarnings <<- FALSE
    }else{
    tryCatch({
        length(hasWarnings)
      },
      error =  function(err){
        wregWarnings <<- list()
        hasWarnings <<- FALSE
        return(hasWarnings)
      })
    }
  }
  
  if(mode =="initialize"){
    init()
  }
  
  if(mode=="add"){
    init()
    wregWarnings <<- c(wregWarnings, w)
    hasWarnings <<- TRUE
  }
  
  if(mode =="get"){
    init()
    if (length(wregWarnings) > 0){
      
      tempWarnings <- wregWarnings
      wregWarnings <<- list()
      hasWarnings <<- FALSE
      
      return(paste("\n",as.character(c("Warning: ",tempWarnings))))
    }
    else{
      return("")
    }
  }
  
  if(mode == "check"){
    init()
    return(hasWarnings)
  }
  
  if(mode == "clear"){
    init(TRUE)
    return(hasWarnings)
  }
  
 
}
