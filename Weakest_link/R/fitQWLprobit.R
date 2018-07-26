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
  require(survival)
  predictor = onedimPredictor(delta, p1, p2)
  if(identical(endpoint, 'ySurv') ) 
    endpoint = with(theData, Surv(time, cens) )
  if(class(endpoint) == 'Surv') {
    result = coxph(endpoint ~ predictor, data=theData)
    #print(result)
    theAIC = 2 - 2*diff(result$loglik)
  }
  else {
    if(class(endpoint)=='character')
      endpoint = theData[[endpoint]]
    if(all(endpoint %in% c(0,1,NA) ) ) 
      fam = binomial
    else 
      fam = normal
    result  = glm(endpoint ~ predictor, 
                  family=fam, data=theData)
    theAIC = result$aic
  }
  if(plotPoints) {
    tryResult = try(points(delta, theAIC, col='blue', pch='X') )
    if(class(tryResult)=='try-error')
      warning('Cannot plot points; probably no active plot.')
  }
  return(list(result=result, theAIC=theAIC))
}

cdf = function(x, FhatStyle=c('normal', 'ecdf')[1]){
  if(FhatStyle == 'normal') 
    return(pnorm(x, mean(x), sd(x)) )
  if(FhatStyle == 'ecdf') 
    return(ecdf(x)(x) )
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
                        covariates = 1, 
                        delta, 
                        dir1 = TRUE, dir2 = TRUE, 
                        testMe = FALSE, 
                        plottheData = TRUE, 
                        FhatStyle=c('normal', 'ecdf')[1],
                        ...) {
  if(testMe)
    theData = WLContinuousdata(...)
  if(length(x1)==1 & is.character(x1)) {
    x1name = x1
    x1 = theData[[x1]]
  }
  else x1name = 'x1'
  if(length(x2)==1 & is.character(x2)) {
    x2name = x2
    x2 = theData[[x2]]
  }
  else x2name = 'x2'
  if(length(endpoint)==1 & is.character(endpoint)) {
    yname = endpoint
    endpoint = theData[[endpoint]]
  }
  else yname = 'y'
  # if(class(endpoint) == 'Surv')
  #   y = theData[[endpoint]]
  if(identical(endpoint, 'ySurv') )
    endpoint = Surv(theData$time, theData$cens)
  #assign('ySurv', attr(theData, 'ySurv'), pos=1, immediate=TRUE)
  Fhat1 = cdf(x1, FhatStyle) * ifelse(dir1, 1, -1)
  Fhat2 = cdf(x2, FhatStyle) * ifelse(dir1, 1, -1)
  if(length(delta) == 1)
    result = fitWithFixedDelta(delta, theData, Fhat1, Fhat2, endpoint)
  else {
      if(missing(delta))
        deltaInterval = c(-3,3)
      else 
        deltaInterval = delta
      optimizeMe = function(delta)
        fitWithFixedDelta(delta, theData, Fhat1, Fhat2, endpoint)$theAIC 
      optresult = optimize(optimizeMe,
                        interval = deltaInterval, 
                        tol = 1e-3)
      optdelta = optresult$minimum #  minimizer!
      optAIC = optresult$objective 
      result = fitWithFixedDelta(optdelta, theData, Fhat1, Fhat2, endpoint)
      attr(result, 'optAIC') = optAIC
      attr(result, 'optdelta') = optdelta
      delta = optdelta   # for plotting
  }  
  
  if(plottheData) {
    #plot(Fhat2, phi2)
    #title(paste('delta = ', signif(digits=3,delta) ) )
    #print(endpoint)
    if(class(endpoint)=='Surv')
      colorChoice =  1+endpoint[ , 'status']
    else
      colorChoice = 1 + (endpoint ) 
    #colorChoice = 1 + (endpoint > median(endpoint)) 
    #print(table(colorChoice))
    plot(x1, x2, pch=c('0','1')[colorChoice], 
         col=c('red','green')[colorChoice],
         xlab = x1name, ylab = x2name)
    drawCOU(x1, x2, delta, FhatStyle) 
    title(paste('delta = ', signif(delta, digits=2)) )
  }
  theframes = sys.frames()
  numframes = length( theframes)
  fr = sys.frame(numframes)
  attr(result, 'frame') = fr
  return(result )
}

#' drawCOU
#' 
#' Draws the COU for a member of the delta family.
#' @details {
#' COU is the set {(x1,x2): phi1(x1) = phi2(x2)}.
#' Equivalently, {(x1,x2): x2 = phi2inv ( phi1(x1) }.
#' 
#' For the quantile stitch weakest link model (QWL),
#' a one-dim predictor is constructed by postulating the locus
#' phi2inv o phi1. 
#' 
#' Quantile stitching says that F2(x2) = F1(x1) on the COU.
#' 
#' So F2(phi2inv ( phi1(F1inv(P)) = P. 
#' 
#'     F2 o phi2inv o phi1 o F1inv = identity.
#'     
#'      phi2inv o phi1 = F2inv o F1. This maps x1 to x2 on the COU.
#'      
#'      x2 = COU(x1) = F2inv ( F1(x1) ).
#'      
#' This COU is in a family of possible COU's indexed by delta:
#' 
#'      x2 = COU(x1) = F2inv ( 1 - H( Hinv(1 - F1(x1)) + delta) )
#'      
#' where H: R-> [0,1]  is symmetric: H(-z) = 1 - H(z)
#' 
#' such as the inverse logit or pnorm.
#' 
#'      x1 = COUinv(x2) = F1inv ( 1 - H( Hinv(1 - F2(x2)) - delta) )
#' }
drawCOU = function(x1, x2, delta, FhatStyle, ...) {
  Fhat1 = cdf(x1, FhatStyle)
  matching_P2 = 1 - H(Hinv(1 - Fhat1) + delta)
  if(FhatStyle == 'normal')
    matching_x2 = qnorm(matching_P2, mean=mean(x2), sd=sd(x2))
  if(FhatStyle == 'ecdf') 
    matching_x2 = quantile(x = x2, probs = matching_P2)
  lines(x1[order(x1)], matching_x2[order(x1)],
        ...)
}