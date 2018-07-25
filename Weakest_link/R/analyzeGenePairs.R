#### analyzeGenePairs: Try all 45 choose 2 gene pairs ####


analyzeAPair = function(g12, g1, g2, endpoint, delta='fit',
                        ...) {
  if(!missing(g12)) {
    g12 = unlist(g12)
    g1 = g12[1]
    g2 = g12[2]
  }
  cat('g1=', g1, ' g2=', g2, '\n')
  d12 = mb[c(g1,g2)]
  #data.frame(x1=mb[[g1]], x2=mb[[g2]])
  endpoint = Surv(mb$time, mb$cens)
  if(identical(delta, 'fit') ) {
    result = fitQWLprobit(testMe = FALSE, x1 = g1, x2 = g2,
                 theData = d12, endpoint = endpoint,
                delta = c(-3,3),
                 ...)
  }
  else result = 
    fitQWLprobit(delta = delta, 
                      theData = d12, endpoint=endpoint,
                      p1 = normcdf(d12$x1), 
                      p2 = normcdf(d12$x2),
                      ...
    )
  return(result)
}

if(interactive()) {
  require(survival)
  allPHmodels = lapply(p45, function(g) {
    theGene = g
    summary(coxph(formula= as.formula(
      paste('Surv(time,cens) ~ ', g)), data=mb, model = F, x = F, y=F)
    )$loglik
  })
  allCoxLogLik = sapply(allPHmodels, diff)
  names(allCoxLogLik) = p45
  summary(allCoxLogLik)
  tail(sort(allCoxLogLik))
  
  #### all pairs ####
  genepairs = expand.grid(sort(p45), sort(p45), 
                          stringsAsFactors = FALSE)
  names(genepairs) = c('g1', 'g2')
  genepairs = genepairs[genepairs$g1 < genepairs$g2, ]
  dim(genepairs)
  choose(45, 2)
  
  #### Test one pair ####
  g1=p45[1]; g2=p45[2]
  randomRow = genepairs[sample(nrow(genepairs), 1), ]
  g1=randomRow$g1; g2=randomRow$g2
  resultFitted = analyzeAPair(g1=g1, g2=g2, 
                        endpoint = mb$D7,  # Surv(mb$time, mb$cens) 
                        delta = 'fit',
                        plottheData= TRUE)
  drawCOU(x1 = mb[[g1]], x2 = mb[[g2]], delta = 0, FhatStyle = 'normal',
          col='blue')
  drawCOU(x1 = mb[[g1]], x2 = mb[[g2]], delta = 0, FhatStyle = 'ecdf',
          col='purple')
  
  resultFixed0 = analyzeAPair(g1=g1, g2=g2, 
               endpoint = mb$D7,  # Surv(mb$time, mb$cens) 
               delta = 0,
               plottheData= TRUE)
  resultFixed0$result$coefficients
  summary(resultFixed0$result)[2]
  pchValue = c(' ', '0', '1')[match(mb$D7, c(NA, 0 , 1))]
  colValue = c('white', 'red', 'blue')[match(mb$D7, c(NA, 0 , 1))]
  plot(mb[[g1]], mb[[g2]], 
       pch = pchValue, col = colValue, cex=0.5)
  
  #### Now, all the pairs.
  allResults = apply(genepairs[1:6, ], 1, function(aRow) {
    g1=aRow[1]; g2=aRow[2]
    analyzeAPair(
        g1=g1, g2=g2, 
        endpoint = mb$D7,  # Surv(mb$time, mb$cens) 
        delta = 'fit',
        plottheData= TRUE)
    drawCOU(x1 = mb[[g1]], x2 = mb[[g2]], delta = 0, 
            FhatStyle = 'normal',
            col='blue')
    drawCOU(x1 = mb[[g1]], x2 = mb[[g2]], delta = 0, 
            FhatStyle = 'ecdf',
            col='purple')
  } )
}