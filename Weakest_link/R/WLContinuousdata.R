WLContinuousdata <-
function(n = 200, p0 = .1, p1 = .5, a1 = 1, a2 = 2, b1 = 1, b2 = 2){
  x1 = rnorm(n = n)
  x2 = rnorm(n = n)
  x1x2 = cbind(x1, x2)
  phi1 = a1 + b1 * x1
  phi2 = a2 + b2 * x2
  phi1phi2 = cbind(phi1, phi2)
  minphi = apply(X = phi1phi2, MARGIN = 1, FUN = min)
  yprob = inverselogit(minphi)
  y = rbinom(n = n, size = 1, prob = yprob)
  return(data.frame(x1, x2, y, minphi))
}