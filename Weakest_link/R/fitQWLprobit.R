#### fitQWLprobit:  fit a parametric weakest link model ####

#' normcdf- empirical CDF fit to a vector x
#'
#' @param x A numeric vector.
#' @return A vector of cumulative probabilities.
#' 
normcdf = function(x) pnorm(x, mean(x), sd(x))

#' compare_cdfs- comparing normcdf to ecdf:
#' 
#' compare_cdfs plots the two cdf's for a particular gene p45[1] CCNB1
#' @details Plots the empirical cdf and the normal-fit cdf.

compare_cdfs = function() {
  data(mb)
  p45x1 = mb[[ p45[1] ]]
  plot( ecdf(x = p45x1) )
  lines(x=qnorm(
    Ptemp<-seq(0,1,length=length(p45x1))
  ), y=Ptemp, col='green', lwd=3)
}

#' deltaMap
#' 
#' Map from [0,1] to [0,1]
#' @param delta A parameter for the map. delta=zero gives the identity function.
#' @param p An input value or vector to the map.
#' 
deltaMap = function(delta, p){
  1 - H(Hinv(1-p) -  delta)
}

#' onedimPredictor
#' 
#' Combine two features into a single predictor.
#' 
#' @param delta An offset for the COU
#' @param p1,p2 Probabilities
#' @return The pmin (parallel min) of p1 and deltaMap(delta, p2) 

onedimPredictor = function(delta, 
                           p1 = Fhat1, p2 = Fhat2){
  if(delta == 0) {
    phi2 <<- p2
  } else if(length(delta)==1) {
    phi2 <<- deltaMap(delta, p2)
  } else
    phi2 <<- sapply(delta, deltaMap, p = p2)
  pmin(p1, phi2)
}

#' onedimPredictorThreeWay
#' 
#' Combine THREE features into a single predictor.
#' 
#' @param delta12,delta13 Offsets for quantile mapping, for forming the COU
#' @param p1,p2,p3 Probabilities
#' @return The pmin (parallel min) of p1 and deltaMap(delta, p2) 
#' @examples {
#' pvec = seq(0,1,length=5)
#' p1p2p3 = expand.grid(pvec,pvec,pvec)
#' result = apply(p1p2p3, 1, function(r) {
#'   p1 = r[1]; p2 = r[2]; p3 = r[3]; 
#'   onedimPredictorThreeWay(0.5, 1, p1=p1, p2=p2, p3=p3)
#' } )
#' p1p2p3r = cbind(p1p2p3, result)
#' names(p1p2p3r) = c('p1','p2','p3', 'result')
#' howManyIsWL = apply(p1p2p3r, 1, function(r){
#'          sum(r[4]==r[1:3]) }  )
#' table(howManyIsWL)
#' whichIsWL = apply(p1p2p3r, 1, function(r){
#'     which(r[4]==r[1:3]) [1]}  )
#' pchWL = c('a','b','c')[whichIsWL]
#' pchWL[howManyIsWL > 1] = 'T'
#' pairs(p1p2p3r, pch = pchWL, col=(2:5)[match(pchWL, c('a','b','c', 'T'))])
#' }

onedimPredictorThreeWay = function(delta12, delta13, 
                           p1 = Fhat1, p2 = Fhat2, p3 = Fhat3){
  if(length(delta12) > 1) stop('delta12 length must be 1')
  if(length(delta13) > 1) stop('delta13 length must be 1')
  if(delta12 == 0) {
    phi2 <<- p2
  } else  {
    phi2 <<- deltaMap(delta12, p2)
  }
  if(delta13 == 0) {
    phi3 <<- p3
  } else  {
    phi3 <<- deltaMap(delta13, p3)
  }
  pmin(p1, phi2, phi3)
}

#' fitWithFixedDelta
#' 
#' Given a value of delta, which determines the COU locus,
#' use coxph or glm to fit 
#' the model with just the one predictor defined by the COU.
#' 
#' @param delta The fixed value of delta, an offset for the COU
#' @param theData Data set.
#' @param p1,p2 CDF values for the two predictors
#' @param endpoint Either the name of the target variable, or the string 'ySurv' to indicate a survival outcome.
#' @param plotPoints If TRUE, add evaluated points to a plot of the AIC vs delta
#' @return A list: 
#' @return  result: The model result outcome object.
#' @return AIC: The AIC or other goodness of fit measure.
#'  
fitWithFixedDelta = function(delta, theData, p1, p2, endpoint, plotPoints = FALSE) {
  predictor = onedimPredictor(delta, p1, p2)
  if(endpoint == 'ySurv') {
    require(survival)
    result = coxph(ySurv ~ predictor, data=theData)
    #print(result)
    theAIC = 2 - 2*diff(result$loglik)
  }
  else {
    if(all(theData[[endpoint]] %in% c(0,1,NA) ) ) 
      fam = binomial
    else 
      fam = normal
    result  = glm(y ~ predictor, family=binomial, data=theData)
    theAIC = result$aic
  }
  if(plotPoints) {
    tryResult = try(points(delta, theAIC, col='blue', pch='X') )
    if(class(tryResult)=='try-error')
      warning('Cannot plot points; probably no active plot.')
  }
  return(list(result=result, theAIC=theAIC))
}

#' fitQWLprobit
#' 
#' Fit a QWL (probit quantile-stitched weakest link model)
#' @param delta The fixed value of delta.
#' @param theData Data set.
#' @param plottheData If TRUE, plot the x1 vs x2 data points, with different characters for y values
fitQWLprobit = function(theData,
                        x1='x1', x2='x2',
                        endpoint='y', ## or 'ySurv'
                        delta, 
                        dir1 = TRUE, dir2 = TRUE, 
                        testMe = FALSE, plottheData = TRUE,
                        ...) {
  if(testMe)
    theData = WLContinuousdata(...)
  x1 = theData[[x1]]
  x2 = theData[[x2]]
  if(endpoint != 'ySurv')
    y = theData[[endpoint]]
  if(endpoint == 'ySurv') 
    assign('ySurv', attr(theData, 'ySurv'), pos=1, immediate=TRUE)
  Fhat1 = pnorm(x1, mean(x1), sd(x1)) * ifelse(dir1, 1, -1)
  Fhat2 = pnorm(x2, mean(x2), sd(x2)) * ifelse(dir2, 1, -1)
  if(length(delta) == 1)
    result = fitWithFixedDelta(delta, theData, Fhat1, Fhat2, endpoint)
  else {
    deltaInterval = ifelse(missing(delta),
                           c(-3,3), delta)
    result = optimize(
      function(delta)
        fitWithFixedDelta(delta, theData, Fhat1, Fhat2, endpoint)$theAIC, 
      interval = deltaInterval, 
      tol = 1e-3)
  }  
  
  if(plottheData) {
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
  theframes = sys.frames()
  numframes = length( theframes)
  fr = sys.frame(numframes)
  attr(result, 'frame') = fr
  return(result )
}

