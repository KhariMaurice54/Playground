#### fitQWLprobit:  fit a parametric weakest link model ####

normcdf = function(x) pnorm(x, mean(x), sd(x))
## comparing normcdf to ecdf:
p45x1 = mb[[ p45[1] ]]
plot( ecdf(x = p45x1) )
lines(x=qnorm(
  Ptemp<-seq(0,1,length=length(p45x1))
), y=Ptemp, col='green', lwd=3)

fitDelta = function(delta, data = data) {
  fitQWLprobit(data = data, testMe = FALSE,
               plotData = FALSE, delta = delta 
  )$aic
}

fitQWLprobit = function(data, endpoint='y',
                        delta, 
                        dir1 = TRUE, dir2 = TRUE, 
                        testMe = FALSE, plotData = TRUE,
                        ...) {
  if(testMe)
    data = WLContinuousdata(...)
  x1 = data$x1
  x2 = data$x2
  y = data$y
  Fhat1 = pnorm(x1, mean(x1), sd(x1)) * ifelse(dir1, 1, -1)
  Fhat2 = pnorm(x2, mean(x2), sd(x2)) * ifelse(dir2, 1, -1)
  H = inverselogit
  Hinv= logit
  fitOneDelta = function(delta) {
    if(delta == 0) {
      phi2 <<- Fhat2
    } else if(length(delta)==1) {
      phi2 <<- 1 - H(Hinv(1-Fhat2) -  delta)
    }
    phiMin = pmin(Fhat1, phi2)
    if(endpoint == 'surv') {
      ySurv = attr(data, 'ySurv')
      require(survival)
      result = coxph(ySurv ~ phiMin)
    }
    else {
      if(all(data[[endpoint]] %in% c(0,1,NA) ) ) 
        fam = binomial
      else 
        fam = normal
      result  = glm(y ~ phiMin, family=binomial, data=data)
    }
  }
  if(length(delta) == 1)
    result = fitOneDelta(delta)
  else {
    deltaInterval = ifelse(missing(delta),
                           c(-3,3), delta)
    result = optimize(fitDelta, 
                      interval = deltaInterval, 
                      tol = 1e-3)
  }  
  
  if(plotData) {
    plot(Fhat2, phi2)
    colorChoice = 1 + (y > median(y))
    plot(x1, x2, pch=c('0','1')[colorChoice], 
         col=c('red','green')[colorChoice])
    # COU is where phi1 = phi2, Fhat1 = phi2,
    #  But phi2 = 1 - H(Hinv(1-Fhat2) -  delta),
    # so Fhat2 = 1 - H(Hinv(1 - Fhat1) + delta)
    matching_P2 = 1 - H(Hinv(1 - Fhat1) + delta)
    matching_x2 = qnorm(matching_P2, mean=mean(x2), sd=sd(x2))
    lines(x1[order(x1)], matching_x2[order(x1)] )
  }
  return(summary(result) )
}

fitQWLprobit(testMe=TRUE, plotData = TRUE, delta=1e-15, 
             b1 = 3, b2 = 5)
testData = WLContinuousdata( b1 = 3, b2 = 5)

deltaSeq = seq(-3,3,length=100)
resultSeq = sapply(deltaSeq, fitDelta, data=testData)
plot(deltaSeq, resultSeq, xlab='delta', ylab='AIC')

optResult = optimize(fitDelta, data=testData,
                     interval = c(-3,3), tol = 1e-3)
abline(v=optResult$minimum, h=optResult$objective, col='red')

