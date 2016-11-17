#quick functions to add and remove warnings

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
