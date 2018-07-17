# metabric_explore\
require(survival)
kmfit = survfit(Surv(time, cens) ~ classif, data=mb)
colors = seq(along=levels(mb$classif))
plot(kmfit, xlab='years', ylab='Survival',
     xscale=365.25, col=colors)
str(mb$classif)
legend(x = 'bottomleft', legend = levels(mb$classif), lty=1, col=colors,
       bty='n', lwd=2)