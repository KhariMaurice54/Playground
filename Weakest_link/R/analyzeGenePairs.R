#### analyzeGenePairs: Try all 45 choose 2 gene pairs ####

makeGenePairs = function() {
  #### all pairs ####
  genepairs = expand.grid(sort(p45), sort(p45), 
                          stringsAsFactors = FALSE)
  names(genepairs) = c('g1', 'g2')
  genepairs = genepairs[genepairs$g1 < genepairs$g2, ]
  genepairs = as.matrix(genepairs)
  return(genepairs)
  # dim(genepairs)
  # choose(45, 2)
}

analyzeAPair = function(g12, g1, g2, endpoint, delta='fit',
                        ...) {
  if(!missing(g12)) {
    g12 = unlist(g12)
    g1 = g12[1]
    g2 = g12[2]
  }
  #cat('g1=', g1, ' g2=', g2, '\n')
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

oneGeneResult = function(geneName='CCNB1', plotMe=TRUE) {
  geneData = mb[[geneName]]
  D7aic = glm(mb$D7 ~ geneData, family=binomial)$aic
  theFormula = Surv(time,cens) ~ geneData
  CoxLogLik = diff(summary(coxph(formula= theFormula, data=mb, 
                            model = F, x = F, y=F)
          )$loglik )
  if(plotMe) {
    boxplot(geneData ~ mb$D7, 
          xlab="died by 7 years", ylab=geneName)
    title(paste(geneName, '\nAIC=', signif(digits=5, D7aic) ))
  }
  return(c(D7aic = D7aic, CoxLL = CoxLogLik) )
}

#### how about single genes? ####
analyzeAllGenes = function(plotMe = FALSE){
  allResultsByGene = sapply(p45, oneGeneResult, plotMe=FALSE)
  allD7AICbyGene <<- allResultsByGene['D7aic', ]
  allCoxLLbyGene <<- allResultsByGene['CoxLL', ]
}

showBestGenes = function() {
  par(mfrow=c(1,1))
  plot(allD7AICbyGene, allCoxLLbyGene)
  lowestAIC = order(allD7AICbyGene)[1]
  g1 = p45[lowestAIC]
  text(allD7AICbyGene[lowestAIC], allCoxLLbyGene[lowestAIC], 
       g1,pos = 4)
  secondhighestLL = order(allCoxLLbyGene,decreasing = TRUE)[2]
  g2 = p45[secondhighestLL]
  text(allD7AICbyGene[secondhighestLL], allCoxLLbyGene[secondhighestLL], 
       g2, pos = 4)
  oneGeneResult(g1, plotMe=TRUE)  # UBE2C
  oneGeneResult(g2, plotMe=TRUE)  # PTTG1
  highestAIC = order(allD7AICbyGene, decreasing = T)[1]
  g3 = p45[highestAIC]
  oneGeneResult(g3, plotMe=TRUE)  # KRT14
}

onePairResult = function(aRow, g1, g2, plottheData= TRUE) {
  if(!missing(aRow)) {
    g1=aRow[1]; g2=aRow[2]
  }
  #cat('g1 ', g1, ' g2 ', g2, '\n')
  result = analyzeAPair(
    g1=g1, g2=g2, 
    endpoint = mb$D7,  # Surv(mb$time, mb$cens) 
    delta = 'fit',
    plottheData= plottheData)
  if(plottheData) {
    drawCOU(x1 = mb[[g1]], x2 = mb[[g2]], delta = 0, 
            FhatStyle = 'normal',
            col='blue')
    drawCOU(x1 = mb[[g1]], x2 = mb[[g2]], delta = 0, 
            FhatStyle = 'ecdf',
            col='purple')
  }
  return(result)
} 
analyzeAllGenePairs = function(plottheData = FALSE,
                               makeBoxPlots = FALSE) {  
  #### Now, all the pairs.
  allD7AICbyPair = apply(genepairs, 1, onePairResult, plottheData = FALSE)
  names(allD7AICbyPair) = apply(genepairs, 1, paste, collapse=',')
  sapply(allD7AICbyPair, `[[`, 'theAIC')
}

showBestPairs = function() {
    par(mfrow=c(1,1))
    lowestAIC = order(allD7AICbyPair)[1]
    g1 = genepairs[lowestAIC, 'g1']
    g2 = genepairs[lowestAIC, 'g2']
    text(allD7AICbyGene[lowestAIC], allCoxLLbyGene[lowestAIC], 
         g1,pos = 4)
    secondhighestLL = order(allCoxLLbyGene,decreasing = TRUE)[2]
    g2 = p45[secondhighestLL]
    text(allD7AICbyGene[secondhighestLL], allCoxLLbyGene[secondhighestLL], 
         g2, pos = 4)
    oneGeneResult(g1, plotMe=TRUE)  # UBE2C
    oneGeneResult(g2, plotMe=TRUE)  # PTTG1
    highestAIC = order(allD7AICbyGene, decreasing = T)[1]
}

plotPairAndSeparateBoxplots = function(g12) {
  #  FOXA1 MYBL2
  g1 = strsplit(g12, ',')[[1]][1]
  g2 = strsplit(g12, ',')[[1]][2]
  result = analyzeAPair(g1=g1, g2=g2)
  par(mfrow=c(1,3))
  boxplot(result$result$linear.predictors ~ mb$D7,
          xlab="died by 7 years", ylab='WL predictor')
  theAIC = title(paste(g1, g2, 
                       '\nAIC=', 
                       signif(digits=4, allD7AICbyGene[g1])
                       #,'   \nCoxLL=', signif(digits=4, allCoxLLbyGene[g1]) 
                 ) )
  oneGeneResult (geneName=g1, plotMe=TRUE) 
  oneGeneResult (geneName=g2, plotMe=TRUE) 
  par(mfrow=c(1,1))
}

if(interactive()) {
  analyzeAllGenes(plotMe = FALSE)
  # produces allD7AICbyGene and  allCoxLLbyGene
  showBestGenes()
  # produces a scatter plot of allD7AICbyGene by allCoxLLbyGene
  # and box plots for the 2 best and one worst.
  analyzeAPair(g12=makeGenePairs()[1,])
  allD7AICbyPair <<- analyzeAllGenePairs(
    plottheData=FALSE, makeBoxPlots = FALSE)
  allD7AICbyPair[paste('PTTG1','UBE2C', sep=',')]
  mean(allD7AICbyPair[paste('PTTG1','UBE2C', sep=',')] <= allD7AICbyPair)
  sum(allD7AICbyPair[paste('PTTG1','UBE2C', sep=',')] > allD7AICbyPair)
  head(sort(allD7AICbyPair))
  sapply(names(head(sort(allD7AICbyPair))),
    plotPairAndSeparateBoxplots
  )
}
