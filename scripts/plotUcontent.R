dataFolder = "/home/alphons/Uni/polar/scripts/u-content/"
range = 200:400
betterScaleStart = 200
betterScaleEnd = 300

tempPos = read.csv(paste0(dataFolder, "resultUcontentPositive.csv"), header=FALSE)
tempNeg = read.csv(paste0(dataFolder, "resultUcontentNegative.csv"), header=FALSE)
tempFN = read.csv(paste0(dataFolder, "resultUcontentFN.csv"), header=FALSE)
tempFP = read.csv(paste0(dataFolder, "resultUcontentFP.csv"), header=FALSE)
tempTP = read.csv(paste0(dataFolder, "resultUcontentTP.csv"), header=FALSE)
tempTN = read.csv(paste0(dataFolder, "resultUcontentTN.csv"), header=FALSE)

xPos = seq(0, 490)
yPos = rep(0, 491)
xNeg = seq(0, 490)
yNeg = rep(0, 491)
xFN = seq(0, 490)
yFN = rep(0, 491)
xFP = seq(0, 490)
yFP = rep(0, 491)
xTP = seq(0, 490)
yTP = rep(0, 491)
xTN = seq(0, 490)
yTN = rep(0, 491)

for (i in 1:(length(yTP) - 1)) {
  yPos[i] = mean(tempPos[, i])
  yNeg[i] = mean(tempNeg[, i])
  yFN[i] = mean(tempFN[, i])
  yTP[i] = mean(tempTP[, i])
  yTN[i] = mean(tempTN[, i])
  yFP[i] = mean(tempFP[, i])
}

par(mfrow=c(3,2), oma=c(0,0,4,0))
plot(xTP[range], yTP[range], main="TP data set", type="h", ylim=c(0,5))
axis(side=1, at=seq(betterScaleStart, betterScaleEnd, by=10))
plot(xTN[range], yTN[range], main="TN data set", type="h", ylim=c(0,5))
axis(side=1, at=seq(betterScaleStart, betterScaleEnd, by=10))
plot(xFN[range], yFN[range], main="FN data set", type="h", ylim=c(0,5))
axis(side=1, at=seq(betterScaleStart, betterScaleEnd, by=10))
plot(xFP[range], yFP[range], main="FP data set", type="h", ylim=c(0,5))
axis(side=1, at=seq(betterScaleStart, betterScaleEnd, by=10))
plot(xPos[range], yPos[range], main="Positive data set", type="h", ylim=c(0,5))
axis(side=1, at=seq(betterScaleStart, betterScaleEnd, by=10))
plot(xNeg[range], yNeg[range], main="Negative data set", type="h", ylim=c(0,5))
axis(side=1, at=seq(betterScaleStart, betterScaleEnd, by=10))
mtext("U content distribution", outer = TRUE, cex = 1.5)