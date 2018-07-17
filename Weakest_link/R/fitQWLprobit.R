#### fitQWLprobit:  fit a parametric weakest link model ####
fitQWLprobit = function(data, delta=0, dir1 = TRUE, dir2 = TRUE, 
                  testMe = FALSE, plotMe = TRUE,
                  ...) {
  if(testMe)
    data = WLContinuousdata(...)
  attach(data)
  on.exit(detach())
  Fhat1 = pnorm(x1, mean(x1), sd(x1)) * ifelse(dir1, 1, -1)
  Fhat2 = pnorm(x2, mean(x2), sd(x2)) * ifelse(dir2, 1, -1)
  if(delta == 0) {
    phi2 = Fhat2
  } else {
    H = inverselogit
    Hinv= logit
    phi2 = 1 - H(Hinv(1-Fhat2) -  delta)
  }
  phiMin = pmin(Fhat1, phi2)
  result  = glm(y ~ phiMin)
  if(plotMe) {
    plot(Fhat2, phi2)
    plot(x1, x2, pch=c('0','1')[y+1], 
         col=c('red','green')[y+1])
    # COU is where phi1 = phi2, Fhat1 = phi2,
    #  But phi2 = 1 - H(Hinv(1-Fhat2) -  delta),
    # so Fhat2 = 1 - H(Hinv(1 - Fhat1) + delta)
    matching_P2 = 1 - H(Hinv(1 - Fhat1) + delta)
    matching_x2 = qnorm(matching_P2, mean=mean(x2), sd=sd(x2))
    lines(x1[order(x1)], matching_x2[order(x1)] )
  }
  summary(result)
}
fitQWLprobit(testMe=TRUE, plotMe = TRUE, delta=2, 
       b1 = 3, b2 = 5)
