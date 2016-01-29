tempTP = read.csv("/home/morikie/Documents/Forschungspraktikum/C++/polar/scripts/resultUcontentTP.csv", header=FALSE)
tempTN = read.csv("/home/morikie/Documents/Forschungspraktikum/C++/polar/scripts/resultUcontentTN.csv", header=FALSE)

xTP = seq(0, 491)
yTP = rep(0, 492)
xTN = seq(0, 491)
yTN = rep(0, 492)

for (i in 1:length(yTP)) {
  yTP[i] = mean(tempTP[, i])
  yTN[i] = mean(tempTN[, i])
}

par(mfrow=c(1,2))
plot(xTP, yTP, main="TP data set", type="h")
plot(xTN, yTN, main="TN data set", type="h")
