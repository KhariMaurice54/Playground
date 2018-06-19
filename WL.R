WLbinarydata = function(n = 200, p0 = .1, p1 = .5){
  x1 = rbinom(n = n, size = 1, prob = 1/2)
  x2 = rbinom(n = n, size = 1, prob = 1/2)
  x1x2 = cbind(x1, x2)
  both = apply(X = x1x2, MARGIN = 1, FUN = min)
  yprob = p0 * (both == 0) + p1 * (both == 1)
  y = rbinom(n = n, size = 1, prob = yprob)
  return(data.frame(x1, x2, y, both))
}
temp=WLbinarydata()
pairs(temp)

inverselogit = function(z){
  exp(z)/ (1 + exp(z))
}
WLContinuousdata = function(n = 200, p0 = .1, p1 = .5, a1 = 1, a2 = 2, b1 = 1, b2 = 2){
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
temp=WLContinuousdata(n = 15)
pairs(temp)

Likelydata = function(data = temp, a1 = 1, a2 = 2, b1 = 1, b2 = 2){
  phij1 = inverselogit(a1 + b1 * data$x1)
  phij2 = inverselogit(a2 + b2 * data$x2)
  phij1phij2 = cbind(phij1, phij2)
  minphij = apply(X = phij1phij2, MARGIN = 1, FUN = min)
  return(sum(data$y*log(minphij) + (1 - data$y) * log(1 - minphij)))
  
}

Likelydata(a1 = .25)
a1seq = seq(0,2,length = 1000)
Lvalues = sapply(a1seq, Likelydata, data = temp)
plot(a1seq, Lvalues, type = "l")
