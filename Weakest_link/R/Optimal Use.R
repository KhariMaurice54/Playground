#The purpose of this function is to plot points on a plane of x1 by x2.
dataset = WLContinuousdata(n = 100, p1 = .9, b1 = 5, b2 = 5)
Plotpoints = function(data = dataset){
  plot(data$x1, data$x2, pch = as.character(data$y), col = c("red", "blue")[1 + data$y])
}
Plotpoints()
attach(as.list(attributes(dataset)[["parameters"]]))
abline(a = ((a1 + a2)/ b2), b = (b1/b2))
detach()
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
  return(list(maximum = Likelydata(a1a2b1b2 = starting_vector), 
           a1maximizer = newa1, 
           a2maximizer = newa2,
           b1maximizer = newb1, 
           b2maximizer = newb2))
}

estimateline = ConMaxStep()
abline(a = (estimateline$a1maximizer - estimateline$a2maximizer)/estimateline$b2maximizer,
       b = estimateline$b1maximizer/estimateline$b2maximizer, 
       col = "purple", 
       lty = 2, 
       lwd = 2) 
