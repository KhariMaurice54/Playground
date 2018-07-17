inverselogit <-
function(z){
  exp(z)/ (1 + exp(z))
}
logit = function(p) log(p/(1-p))
