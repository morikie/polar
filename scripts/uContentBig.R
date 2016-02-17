bigUContentTN = read.csv("/home/morikie/Documents/Forschungspraktikum/C++/polar/scripts/uContentBigTN.txt", header=FALSE)
bigUContentTP = read.csv("/home/morikie/Documents/Forschungspraktikum/C++/polar/scripts/uContentBigTP.txt", header=FALSE)

TNdistribution = sort(bigUContentTN[, 1])
TPdistribution = sort(bigUContentTP[, 1])

hist(TNdistribution)
hist(TPdistribution)

#yTN = 1:length(TNdistribution)
#yTP = 1:length(TPdistribution)

#par(mfrow=c(1,2))
#plot(yTN, TNdistribution)
#plot(yTP, TPdistribution)