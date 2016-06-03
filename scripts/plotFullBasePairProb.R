tempPos = read.csv("/home/morikie/Documents/Forschungspraktikum/C++/polar/scripts/utr_foldings/fullBasePairProb.csv", header=FALSE)

numPlots = 10

par(mfrow=c(5,2), oma=c(0,0,4,0))

for (k in sample(length(tempPos[, 1]), numPlots, replace=F)) {
  
  len = 0
  for (i in (length(tempPos[k,] - 1)):1) {
    if (! is.na(tempPos[k, i])) {
      len = i
      break
    }
  }
  
  xPos = seq(0, len - 5)
  yPos = rep(0, len - 4)
  pasPos = rep(0, 5)
  
  for (i in 1:5) {
    pasPos[i] = as.integer(tempPos[k, i])
  }
  for (i in 6:(len + 1)) {
    if (! is.na(tempPos[k, i])) {
      yPos[i] = tempPos[k, i]
    } else {
      yPos[i] = 0.0 
    }
    
  }
  
  plot(xPos, yPos, main=paste("Bpp of row ", toString(k)), type="h", ylim=c(0,1))
  
  for (i in 1:5) {
    if (pasPos[i] != -1) {
      abline(v=pasPos[i], col="red")
    }
  }
}
mtext("Bpp of random 3' UTR sequences", outer = TRUE, cex = 1.5)