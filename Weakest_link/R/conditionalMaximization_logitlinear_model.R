####Two Dimensional Plot####
Two_Dimensional_Plot = function( 
  dataset = WLContinuousdata(),
  gridSize = 200,
  a1seq = seq(0,2,length = gridSize),
  a2seq = seq(0,4,length = gridSize)
){
  parameters = attr(dataset, "parameters")
  a1 = parameters['a1']
  a2 = parameters['a2']
  b1 = parameters['b1']
  b2 = parameters['b2']
  a1a2grid = expand.grid(a1seq, a2seq)
  Lvalues = apply(a1a2grid, 1, function(row){
    Likelydata(data = dataset, a1 = row[1], a2 = row[2])})
  Lvalues = matrix(Lvalues, nrow = length(a1seq))
  contour.default(x = a1seq, y =a2seq, z = Lvalues,
                  xlab='a1', ylab='a2')
  title('Log likelihood surface')
  points(a1, a2, col = "green", pch = "x", cex = 2)
  maximizer = which(max(Lvalues) == Lvalues)
  rmax = row(Lvalues)[maximizer]
  cmax = col(Lvalues)[maximizer]
  points(a1seq[rmax], a2seq[cmax], col = "red", pch = "o", cex = 2)
  legend(x = 'top', col=c('red','green'),
         pch=c('o', 'x'), legend = c('grid max', 'true a1,a2'))
  return(c(maximum = max(Lvalues), a1maximizer = a1seq[rmax], a2maximizer = a2seq[cmax]))
}
if(interactive())
  Two_Dimensional_Plot()

####Finding a more accurate maximizer liklihood estimate by "hem-stiching"(conditional maximization)####

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
# ConMaxResult = ConMax(tol = 1e-9)
# pairs(ConMaxResult)
