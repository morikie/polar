tempTP = read.csv("/home/morikie/Documents/Forschungspraktikum/C++/polar/scripts/resultUcontentTP.csv", header=FALSE)
tempTN = read.csv("/home/morikie/Documents/Forschungspraktikum/C++/polar/scripts/resultUcontentTN.csv", header=FALSE)
tempFN = read.csv("/home/morikie/Documents/Forschungspraktikum/C++/polar/scripts/notFoundPAS.csv", header=FALSE)

xFN = seq(0, 491)
yFN = rep(0, 492)
xTP = seq(0, 491)
yTP = rep(0, 492)
xTN = seq(0, 491)
yTN = rep(0, 492)

for (i in 1:(length(yTP) - 1)) {
  yFN[i] = mean(tempFN[, i])
  yTP[i] = mean(tempTP[, i])
  yTN[i] = mean(tempTN[, i])
}

par(mfrow=c(1,3))
plot(xFN[200:300], yFN[200:300], main="FN data set", type="h", ylim=c(0,4))
axis(side=1, at=seq(200, 300, by=10))
plot(xTP[200:300], yTP[200:300], main="TP data set", type="h", ylim=c(0,4))
axis(side=1, at=seq(200, 300, by=10))
plot(xTN[200:300], yTN[200:300], main="TN data set", type="h", ylim=c(0,4))
axis(side=1, at=seq(200, 300, by=10))



