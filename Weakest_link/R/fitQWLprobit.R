#### fitQWLprobit:  fit a parametric weakest link model ####

normcdf = function(x) pnorm(x, mean(x), sd(x))

compare_cdfs = function() {
  ## comparing normcdf to ecdf:
  data(mb)
  p45x1 = mb[[ p45[1] ]]
  plot( ecdf(x = p45x1) )
  lines(x=qnorm(
    Ptemp<-seq(0,1,length=length(p45x1))
  ), y=Ptemp, col='green', lwd=3)
}

deltaMap = function(delta, p){
  1 - H(Hinv(1-p) -  delta)
}

onedimPredictor = function(delta, 
                           p1 = Fhat1, p2 = Fhat2){
  if(delta == 0) {
    phi2 <<- p2
  } else if(length(delta)==1) {
    phi2 <<- deltaMap(delta, p2)
  }
  pmin(p1, phi2)
}
fitOneDelta = function(delta, p1, p2, endpoint) {
  predictor = onedimPredictor(delta, p1, p2)
  if(endpoint == 'ySurv') {
    require(survival)
    result = coxph(ySurv ~ predictor)
    #print(result)
    theAIC = 2 - 2*diff(result$loglik)
  }
  else {
    if(all(data[[endpoint]] %in% c(0,1,NA) ) ) 
      fam = binomial
    else 
      fam = normal
    result  = glm(y ~ predictor, family=binomial, data=data)
    theAIC = result$aic
  }
  return(list(result=result, theAIC=theAIC))
}

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
  if(length(delta) == 1)
    result = fitOneDelta(delta, Fhat1, Fhat2, endpoint)
  else {
    deltaInterval = ifelse(missing(delta),
                           c(-3,3), delta)
    result = optimize(
      function(delta)
        fitDelta(delta, Fhat1, Fhat2, endpoint)$theAIC, 
      interval = deltaInterval, 
      tol = 1e-3)
  }  
  
  if(plotData) {
    plot(Fhat2, phi2)
    print(endpoint)
    if(endpoint=='ySurv')
      colorChoice =  1+ySurv[ , 'status']
    else
      colorChoice = 1 + (y > median(y)) 
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
  attr(result, 'frame') = sys.frame(1)
  return(result )
}

fitQWLprobitTests = function() {
  fitQWLprobit(testMe=TRUE, plotData = FALSE, delta=1e-15, 
               b1 = 3, b2 = 5)
  testData = WLContinuousdata( b1 = 3, b2 = 5)
  deltaSeq = seq(-3,3,length=100)
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
                       interval = c(-3,3), tol = 1e-7)
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
  
  #### evaluation ####
  #### prediction from a fitQWLprobit model
  # We calculate the onedimPredictor
  testResult = fitQWLprobit(testMe=TRUE, plotData = FALSE, delta=1e-15, 
                            b1 = 3, b2 = 5)
  theFrame = attr(testResult, 'frame')
  ls(env=theFrame)
  #attach(theFrame)
  with(theFrame, {
    H=get('H', env=theFrame)
    onedimPredictor(delta = delta)
  }
  )
}
