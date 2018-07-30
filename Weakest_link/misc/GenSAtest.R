library(GenSA)
sv = unlist(starting_vector)
SA = GenSA(par = sv, fn = NegaLike, lower = sv - 3, upper = sv + 3, control = list(maxit = 100))
drawCOUline(a1 = SA$par["a1"], a2 = SA$par["a2"], b1 = SA$par['b1'], b2 = SA$par["b2"])
