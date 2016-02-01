tempTP = read.csv("/home/morikie/Documents/Forschungspraktikum/C++/polar/scripts/resultUcontentTP.csv", header=FALSE)
tempTN = read.csv("/home/morikie/Documents/Forschungspraktikum/C++/polar/scripts/resultUcontentTN.csv", header=FALSE)

xTP = seq(0, 490)
yTP = rep(0, 491)
xTN = seq(0, 490)
yTN = rep(0, 491)

for (i in 1:length(yTP)) {
  yTP[i] = mean(tempTP[, i])
  yTN[i] = mean(tempTN[, i])
}

par(mfrow=c(1,2))
plot(xTP[240:260], yTP[240:260], main="TP data set", type="h")
plot(xTN[240:260], yTN[240:260], main="TN data set", type="h")
