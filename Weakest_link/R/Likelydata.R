Likelydata <-
function(data = temp, a1 = 1, a2 = 2, b1 = 1, b2 = 2){
  phij1 = inverselogit(a1 + b1 * data$x1)
  phij2 = inverselogit(a2 + b2 * data$x2)
  phij1phij2 = cbind(phij1, phij2)
  minphij = apply(X = phij1phij2, MARGIN = 1, FUN = min)
  return(sum(data$y*log(minphij) + (1 - data$y) * log(1 - minphij)))
  
}
