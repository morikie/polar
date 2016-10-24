dataFolder = "/home/alphons/Uni/polar/scripts/"
tempPos = read.csv(paste0(dataFolder, "utr_foldings/basePairProbAtPAS.csv"), header=FALSE)

xPos = seq(0, 99)
yPos = rep(0, 100)

for (i in 1:(length(yPos) - 1)) {
  yPos[i] = mean(tempPos[, i])
}

plot(xPos, yPos, main="Base pair probability data set", type="h", ylim=c(0,0.5))
axis(side=1, at=seq(0, 100, by=5))
abline(v=74, col="red")
abline(v=81, col="red")
