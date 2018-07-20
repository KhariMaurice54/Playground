.onAttach = function(libname, pkgname) {
  data(mb)
  attr(mb, 'ySurv') <<- with(data=mb, survival::Surv(time, cens))
  # double-headed <<- is required
  data(p45)
  data(PAM50genes)
}