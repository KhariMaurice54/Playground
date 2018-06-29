require(WeakestLink)
temp= WLbinarydata()
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
B1_Varies = function(b1){
  Likelydata(b1 = b1, 
             a1 = starting_vector$a1, 
             a2 = starting_vector$a2, 
             b2 = starting_vector$b2)
}
B2_Varies = function(b2){
  Likelydata(a1 = starting_vector$a1, 
             b1 = starting_vector$b1, 
             b2 = b2, 
             a2 = starting_vector$a2)
}
ConMaxStep = function(){
  newa1 = optimize(f = A1_Varies, interval = c(0,2), maximum = TRUE)$maximum
  starting_vector$a1 <<- newa1
  newa2 = optimize(f = A2_Varies, interval = c(0,2), maximum = TRUE)$maximum
  starting_vector$a2 <<- newa2
  newb1 = optimize(f = B1_Varies, interval = c(0,2), maximum = TRUE)$maximum
  starting_vector$b1 <<- newb1
  newb2 = optimize(f = B2_Varies, interval = c(0,2), maximum = TRUE)$maximum
  starting_vector$b2 <<- newb2
  return(c(maximum = Likelydata(a1a2b1b2 = starting_vector), 
           a1maximizer = newa1, 
           a2maximizer = newa2,
           b1maximizer = newb1, 
           b2maximizer = newb2))
}
ConMax = function(tol = 1e-7, 
                  starting_vector = list(a1 = 1, a2 = 3, b1 = 0, b2 = 1)
){
  starting_vector <<- starting_vector
  shouldstop = FALSE
  Results = NULL
  while(shouldstop == FALSE){
    Temp = ConMaxStep()
    if(is.null(Results)){
      Results = t(as.data.frame(Temp))
    }
    else{
      Results = rbind(Results, Temp)
      Delta = Results[nrow(Results), "maximum"] -
        Results[nrow(Results) - 1, "maximum"]
      if(Delta < tol){
        shouldstop = TRUE
      }
    }
  }
  return(Results)
}
ConMaxResult = ConMax(tol = 1e-9)
pairs(ConMaxResult)
