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

#Two Dimensional Plot
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
}
Two_Dimensional()


