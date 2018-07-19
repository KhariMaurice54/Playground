#' Draw a curve of optimal use from a parametric model.
#' 
#' \code{drawCOUline} Draws a 'curve of optimal use' (COU) line on the current plot.
#' #' 
#' @param a1 intercept for inverselogit(phi1(x1))
#' @param a2 intercept for inverselogit(phi2(x2))
#' @param b1 slope for inverselogit(phi1(x1))
#' @param b2 slope for inverselogit(phi2(x2))
#' @param ... Parameters passed to \code{abline()}
#' 
#' @details 
#' solves the equation \code{a1+b1*x1 = a2+b2*x2} for x2.
#' 
drawCOUline = function(a1, a2, b1, b2, ...){
  abline(a = ((a1 + a2)/ b2), b = (b1/b2), ...)
}
