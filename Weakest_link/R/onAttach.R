.onAttach = function(libname, pkgname) {
  data(p45)
  genepairs <<- makeGenePairs()
  data(PAM50genes)
  data(mb)
  mb$D7 <<- (mb$time < 7*365.25) & (mb$cens == 1)
  mb$D7[(mb$time < 7*365.25) & (mb$cens == 0)] <<- NA
  attr(mb, 'ySurv') <<- with(data=mb, survival::Surv(time, cens))
  D7label = 'Died by 7 years'
  # double-headed <<- is required
}
