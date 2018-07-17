# metabric_explore\
require(survival)
kmfit = survfit(Surv(time, cens) ~ classif, data=mb)
colors = seq(along=levels(mb$classif))
plot(kmfit, xlab='years', ylab='Survival',
     xscale=365.25, col=colors)
str(mb$classif)
legend(x = 'bottomleft', legend = levels(mb$classif), lty=1, col=colors,
       bty='n', lwd=2)

p45test = sapply(p45, function(g) 
  {
  #cat('=== ', g, ' ===\n')
  g = mb[[g]]
  result = summary(coxph(Surv(time, cens) ~ g, data=mb) ) 
  c(logtest=result$logtest['pvalue'], sctest=result$sctest['pvalue'])
})
p45plog = log10(p45test['logtest.pvalue', ])
plot.ecdf(p45plog)
points(sort(p45plog), (1:45)/45, col='blue', pch='x')
text(sort(p45plog) [1], ( (1:45)/45) [1], 
     p45[order(p45plog)] [1], pos = 3)

plot(survfit(Surv(time, cens) ~ (UBE2C > median(UBE2C)), 
             data=mb) )
with(data=mb,
     table(classif, (UBE2C > median(UBE2C)) 
     ) 
)
