#' Utility function to format the regression equation
#' 
#' @description
#' The \code{regEquationFormat} function formats 
#' the regression equation for WREG-R.
#' 
#' @param X The X/Y variable
#' @param C1 (optional) Constant 1 for transform
#' @param C2 (optional) Constant 2 for transform
#' @param C3 (optional) Constant 3 for transform
#' @param C4 (optional) Constant 4 for transform
#' @param mode Type of transform on X/Y variable
#' @param var Whether or not method applies to the X or Y variable
#' 
#' @details
#' This functions formats the regression equation based on input for WREG output.
#' 
#' @return Nothing, strings are added to list xEq
#' 
#'@examples
#'regEquationFormat("P", mode="none", var="X")
#'
#'@export

regEquationFormat <- function(X, C1 = NULL, C2 = NULL, C3 = NULL, C4 = NULL, mode, var){
  if(mode == "none"){
    if(var == "X"){
      addXEquation(X)
    }else{
      yEq <<- sprintf("%s = ", X)
    }
   
  }
  else{
    part1 <- ifelse(mode=="log10","log_{10}",mode)
    bracket1 <- "{("
    c1 <- ifelse(C1 == 1,"",sprintf("%.2f * ", C1))
    c2 <- ifelse(C2 == 0 | C2 == 1,X, sprintf("%s^{%.2f}", X, C2))
    c3 <- ifelse(C3 == 0,"",sprintf(" + %.2f", C3))
    bracket2 <- ")"
    c4 <- ifelse(C4 == 0 | C4 == 1,"",sprintf("^{%.2f}", C4))
    bracket3 <- "}"
    
    if (var == "X"){
      addXEquation(paste(c(part1, bracket1, c1, c2, c3, bracket2, c4, bracket3), collapse=""))
    }else{
      yEq <<- sprintf("%s = ", paste(c(part1, bracket1, c1, c2, c3, bracket2, c4, bracket3), collapse=""))
    }
    
  }
}

addXEquation = function(equation){
  if(is.null(xEq)){
    xEq <<- equation
  }else{
    xEq <<- c(xEq, equation)
  }
}