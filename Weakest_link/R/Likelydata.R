#' @param data temporary dataset the function is working with
#' @param a1 intercept for inverselogit(phi2(x2))
#' @param a2 slope for inverselogit(phi1(x1))
#' @param b1 slope for inverselogit(phi2(x2))
#' @param b2 slope for inverselogit(phi2(x2))
#' @param ... Parameters passed to \code{abline()}
Likelydata <-
function(data = temp, a1 = 1, a2 = 2, b1 = 1, b2 = 2, a1a2b1b2){
  if(!missing(a1a2b1b2)){
    a1a2b1b2 = unlist(a1a2b1b2)
    a1 = a1a2b1b2['a1']
    a2 = a1a2b1b2['a2']
    b1 = a1a2b1b2['b1']
    b2 = a1a2b1b2['b2']
  }
  phij1 = inverselogit(a1 + b1 * data$x1)
  phij2 = inverselogit(a2 + b2 * data$x2)
  phij1phij2 = cbind(phij1, phij2)
  minphij = apply(X = phij1phij2, MARGIN = 1, FUN = min)
  return(sum(data$y*log(minphij) + (1 - data$y) * log(1 - minphij)))
  
}

NegaLike = function(par){
  ### For use as input to GenSA, which only does minimization.
  (-1)*Likelydata(a1a2b1b2 = as.list(par))
  
}
