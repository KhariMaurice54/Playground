#### fitQWLprobit:  fit a parametric weakest link model ####

normcdf = function(x) pnorm(x, mean(x), sd(x))
## comparing normcdf to ecdf:
p45x1 = mb[[ p45[1] ]]
plot( ecdf(x = p45x1) )
lines(x=qnorm(
  Ptemp<-seq(0,1,length=length(p45x1))
), y=Ptemp, col='green', lwd=3)

fitDelta = function(delta, plotPoints = FALSE,...) {
  result = fitQWLprobit(testMe = FALSE,
               plotData = FALSE, delta = delta,
               ...
  )$theAIC
  if(plotPoints)
    points(delta, result, col='blue', pch='X')
  return(result)
}

fitQWLprobit = function(data,
                        x1='x1', x2='x2',
                        endpoint='y', ## or 'ySurv'
                        delta, 
                        dir1 = TRUE, dir2 = TRUE, 
                        testMe = FALSE, plotData = TRUE,
                        ...) {
  if(testMe)
    data = WLContinuousdata(...)
  x1 = data[[x1]]
  x2 = data[[x2]]
  if(endpoint != 'ySurv')
    y = data[[endpoint]]
  if(endpoint == 'ySurv') 
    assign('ySurv', attr(data, 'ySurv'), pos=1, immediate=TRUE)
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
    if(endpoint == 'ySurv') {
      require(survival)
      result = coxph(ySurv ~ phiMin)
      #print(result)
      theAIC = 2 - 2*diff(result$loglik)
    }
    else {
      if(all(data[[endpoint]] %in% c(0,1,NA) ) ) 
        fam = binomial
      else 
        fam = normal
      result  = glm(y ~ phiMin, family=binomial, data=data)
      theAIC = result$aic
    }
    return(list(result=result, theAIC=theAIC))
  }
  if(length(delta) == 1)
    result = fitOneDelta(delta)
  else {
    deltaInterval = ifelse(missing(delta),
                           c(-3,3), delta)
    result = optimize(
      function(delta)
        fitDelta(delta)$theAIC, 
      interval = deltaInterval, 
      tol = 1e-3)
  }  
  
  if(plotData) {
    plot(Fhat2, phi2)
    print(endpoint)
    colorChoice = if(endpoint=='ySurv') 
                         1+ySurv[ , 'status'] else
                         1 + (y > median(y)) 
    print(table(colorChoice))
    plot(x1, x2, pch=c('0','1')[colorChoice], 
         col=c('red','green')[colorChoice])
    # COU is where phi1 = phi2, Fhat1 = phi2,
    #  But phi2 = 1 - H(Hinv(1-Fhat2) -  delta),
    # so Fhat2 = 1 - H(Hinv(1 - Fhat1) + delta)
    matching_P2 = 1 - H(Hinv(1 - Fhat1) + delta)
    matching_x2 = qnorm(matching_P2, mean=mean(x2), sd=sd(x2))
    lines(x1[order(x1)], matching_x2[order(x1)] )
  }
  
  return(result )
}

fitQWLprobit(testMe=TRUE, plotData = TRUE, delta=1e-15, 
             b1 = 3, b2 = 5)
testData = WLContinuousdata( b1 = 3, b2 = 5)

deltaSeq = seq(0,3,length=100)
resultSeq = sapply(deltaSeq, fitDelta, data=testData)
plot(deltaSeq, resultSeq, xlab='delta', ylab='AIC')

optResult = optimize(fitDelta, data=testData,
                     interval = c(-3,3), tol = 1e-3)
abline(v=optResult$minimum, h=optResult$objective, col='red')

##### Now on real survival data ####
fitQWLprobit(data = mb, delta=0,
             x1=names(sort(p45plog))[1],
             x2=names(sort(p45plog))[2],
             endpoint='ySurv')
deltaSeq = seq(0,3,length=500)
resultSeq = sapply(deltaSeq, fitDelta, 
                   data=mb, 
                   x1=names(sort(p45plog))[1],
                   x2=names(sort(p45plog))[2],
                   endpoint='ySurv' )
plot(deltaSeq, resultSeq, xlab='delta', ylab='AIC')

optResult = optimize(fitDelta, data=mb, 
                     x1=names(sort(p45plog))[1],
                     x2=names(sort(p45plog))[2],
                     endpoint='ySurv',
                     interval = c(0,3), tol = 1e-7)
abline(v=optResult$minimum, h=optResult$objective, col='red')

install.packages('GenSA')
help(p='GenSA')
require(GenSA)
system.time(saResult<<-GenSA(par=0, lower=-1, upper= 2,
      control=list(maxit=1000),
      fitDelta, plotPoints = TRUE, data=mb, 
      x1=names(sort(p45plog))[1],
      x2=names(sort(p45plog))[2],
      endpoint='ySurv'))
str(saResult)
plot(saResult$trace.mat[,1], saResult$trace.mat[,2], log='y')
plot(saResult$trace.mat[,1], saResult$trace.mat[,3],
     ylim=c(saResult$value, -92.5))
plot(saResult$trace.mat[,1], saResult$trace.mat[,4],
     ylim=c(saResult$value, -92.7))
plot(saResult$trace.mat[,3], saResult$trace.mat[,4])

