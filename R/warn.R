#quick functions to add and remove warnings

warn <- function(mode, w=NULL){

  init <- function(){
    tryCatch({
      length(hasWarnings)
    },
    error =  function(err){
      wregWarnings <<- list()
      hasWarnings <<- FALSE
      return(hasWarnings)
    })
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
  
  
 
}
