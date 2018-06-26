require(WeakestLink)
temp=WLbinarydata()
pairs(temp)

temp=WLContinuousdata(n = 150)
pairs(temp)


Likelydata()
a1seq = seq(0,2,length = 100)
a2seq = seq(0,4,length = 100)
#One dimensional plot
Lvalues = sapply(a1seq, Likelydata, data = temp)
plot(a1seq, Lvalues, type = "l")
abline(v = a1seq[which.max(Lvalues)])
abline(h = max(Lvalues))

####Two Dimensional Plot####
Two_Dimensional = function(){
  a1a2grid = expand.grid(a1seq, a2seq)
  Lvalues = apply(a1a2grid, 1, function(row){Likelydata(a1 = row[1], a2 = row[2])})
  Lvalues = matrix(Lvalues, nrow = 100)
  contour.default(x = a1seq, y =a2seq, z = Lvalues)
  args(WLContinuousdata)
  points(1, 2, col = "green", pch = "x", cex = 2)
  maximizer = which(max(Lvalues) == Lvalues)
  rmax = row(Lvalues)[maximizer]
  cmax = col(Lvalues)[maximizer]
  points(a1seq[rmax], a2seq[cmax], col = "red", pch = "o", cex = 4)
  return(c(maximum = max(Lvalues), a1maximizer = a1seq[rmax], a2maximizer = a2seq[cmax]))
}
Two_Dimensional()

####Finding a more accurate maximizer liklihood estimate by "hem-stiching"(conditional maximization)####
Four_Dimensional = function(starting_vector = list(a1 = 1, a2 = 1, b1 = 1, b2 = 1)){
  optimize()
}

####A1 Variation Function####
A1_Varies = function(a1){
  Likelydata(a1 = a1, 
             b1 = starting_vector$b1, 
             a2 = starting_vector$a2, 
             b2 = starting_vector$b2)
}
A2_Varies = function(a2){
  Likelydata(a1 = starting_vector$a1, 
             b1 = starting_vector$b1, 
             a2 = a2, 
             b2 = starting_vector$b2)
}
starting_vector = list(a1 = 0, a2 = 0, b1 = 1, b2 = 2)
ConMaxResult = sapply(1:10, function(n){
  newa1 = optimize(f = A1_Varies, interval = c(0,2), maximum = TRUE)$maximum
  newa2 = optimize(f = A2_Varies, interval = c(0,2), maximum = TRUE)$maximum
  starting_vector$a1 <<- newa1
  starting_vector$a2 <<- newa2
  points(newa1, newa2, col = "green", pch = "o", cex = 4)
  return(c(maximum = Likelydata(starting_vector), a1maximizer = newa1, a2maximizer = newa2))
}
)
Likelydata(a1 = 0.7351119, a2 =0.9536362, b1 = 1, b2 = 2)

