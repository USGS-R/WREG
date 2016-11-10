#quick functions to validate inputs
#returns FALSE for warnings if the input is valid and TRUE if the if the input is not valid

wregValidation <- function(input, mode, compare=NULL, message=NULL, warnFlag = FALSE){

  invalid <- FALSE
  
  error <- function(default){
    print(warnFlag)
    if(warnFlag){
      #print for command line
      print(paste0("Warning ",ifelse(
        is.null(message),default,message)))
      
      #add warning
      warn("add", ifelse(
        is.null(message),default,message))
      invalid <- TRUE
    }else{
      #halt program (to be caught or stopped in command line)
      stop(ifelse(
        is.null(message),default,message))
    }
  }
  
  numeric <- function(x) {
    if(is.list(x)){
      if(!Reduce('&',lapply(x, is.numeric))){
        error("Input needs to be numeric")
      }
    }else{
      if(!is.numeric(x)){
        error("Input needs to be numeric")
      }
    }
  }
  
  greaterThan <- function(x){
    if(is.list(x)){
      if(!Reduce('&',lapply(x, function(y) {return(y > compare)}))){
        error(paste0("Input needs to be greater than ", compare))
      } 
    }
    else{
      if(!(x > compare)){
        error(paste0("Input needs to be greater than ", compare))
      }
    }
  }
  
  greaterThanEq <- function(x){
    if(is.list(x)){
      if(!Reduce('&',lapply(x, function(y) {return(y >= compare)}))){
        error(paste0("Input needs to be greater than equal to", compare))
      } 
    }
    else{
      if(!(x >= compare)){
        error(paste0("Input needs to be greater than equal to", compare))
      }
    }
  }
  
  lessThan <- function(x){
    if(is.list(x)){
      if(!Reduce('&',lapply(x, function(y) {return(y < compare)}))){
        error(paste0("Input needs to be less than ", compare))
      } 
    }
    else{
      if(!(x < compare)){
        error(paste0("Input needs to be less than ", compare))
      }
    }
  }
  
  lessThanEq <- function(x){
    if(is.list(x)){
      if(!Reduce('&',lapply(x, function(y) {return(y <= compare)}))){
        error(paste0("Input needs to be less than ", compare))
      } 
    }
    else{
      if(!(x <= compare)){
        error(paste0("Input needs to be less than or equal to ", compare))
      }
    }
  }
  
  eq <- function(x){
    if(is.list(x)){
      if(!Reduce('&',lapply(x, function(y) {return(y == compare)}))){
        error(paste0("Input needs to be equal to ", compare))
      } 
    }
    else{
      if(!(x == compare)){
        error(paste0("Input needs to be equal to ", compare))
      }
    }
  }
  
  notEq <- function(x){
    if(is.list(x)){
      if(!Reduce('&',lapply(x, function(y) {return(y != compare)}))){
        error(paste0("Input cannot be equal to ", compare))
      } 
    }
    else{
      if(!(x != compare)){
        error(paste0("Input cannot be equal to ", compare))
      }
    }
  }
  
  if(mode == "greaterThan"){
    numeric(input)
    greaterThan(input)
    return(invalid)
    
  }
  
  if(mode == "greaterThanEq"){
    numeric(input)
    greaterThanEq(input)
    return(invalid)
  }
  
  if(mode == "lessThan"){
    numeric(input)
    lessThan(input)
    return(invalid)
  }
  
  if(mode == "lessThanEq"){
    numeric(input)
    lessThanEq(input)
    return(invalid)
  }
  
  if(mode == "numeric"){
    numeric(input)
    return(invalid)
  }
  
  if(mode == "eq"){
    eq(input)
    return(invalid)
  }
  
  if(mode == "notEq"){
    notEq(input)
    return(invalid)
  }
  
}