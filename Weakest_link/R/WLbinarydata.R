WLbinarydata <-
function(n = 200, p0 = .1, p1 = .5){
  x1 = rbinom(n = n, size = 1, prob = 1/2)
  x2 = rbinom(n = n, size = 1, prob = 1/2)
  x1x2 = cbind(x1, x2)
  both = apply(X = x1x2, MARGIN = 1, FUN = min)
  yprob = p0 * (both == 0) + p1 * (both == 1)
  y = rbinom(n = n, size = 1, prob = yprob)
  return(data.frame(x1, x2, y, both))
}
