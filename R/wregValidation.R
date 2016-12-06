#' Validation Scheme for WREG
#' 
#' @description
#' The \code{wregValidation} function validates input in WREG-R.
#' 
#' @param input A value to validate.
#' @param mode Method of validation
#' @param compare (optional) Value to compare against if mode requires comparison
#' @param message (optional) Custom message for warning/error
#' @param warnFlag (optional) Whether or not invalid input should throw a warning/error
#' 
#' @details
#' This functions streamlines validation for the programmer so less lines of code are written
#' 
#' @return TRUE or FALSE depending on validity of input or throws an error
#' 
#'@examples
#'invalid <- wregValidation(0, "eq", 0)
#'
#'@export

wregValidation <- function(input, mode, compare=NULL, message=NULL, warnFlag = FALSE){

  invalid <<- FALSE
  
  error <- function(default){
    if(warnFlag){
      #print for command line
      print(paste0("Warning ",ifelse(
        is.null(message),default,message)))
      
      #add warning
      warn("add", ifelse(
        is.null(message),default,message))
      invalid <<- TRUE
    }else{
      #halt program (to be caught or stopped in command line)
      stop(ifelse(
        is.null(message),default,message))
    }
  }
  
  numeric <- function(x) {
    if(length(x) > 1){
      if(!Reduce('&',lapply(x, is.numeric))){
        error("Input needs to be numeric")
      }
    }else{
      if(!is.numeric(x)){
        error("Input needs to be numeric")
      }
    }
  }
  
  infinite <- function(x) {
    if(length(x) > 1){
      if(Reduce('|',lapply(x, is.infinite))){
        error("Input cannot be infinite")
      }
    }else{
      if(is.infinite(x)){
        error("Input cannot be infinite")
      }
    }
  }
  
  greaterThan <- function(x){
    if(length(x) > 1){
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
    if(length(x) > 1){
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
    if(length(x) > 1){
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
    if(length(x) > 1){
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
    if(length(x) > 1){
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
    if(length(x) > 1){
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
  
  if(mode == "infinite"){
    infinite(input)
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