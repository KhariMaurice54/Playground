inverselogit <-
function(z){
  exp(z)/ (1 + exp(z))
}
