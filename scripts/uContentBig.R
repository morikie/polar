bigUContentTN = read.csv("/home/alphons/Uni/Forschungspraktikum/polar/scripts/u-content/uContentBigTN.txt", header=FALSE)
bigUContentTP = read.csv("/home/alphons/Uni/Forschungspraktikum/polar/scripts/u-content/uContentBigTP.txt", header=FALSE)
smallUContentTN = read.csv("/home/alphons/Uni/Forschungspraktikum/polar/scripts/u-content/uContentSmallTN.txt", header=FALSE)
smallUContentTP = read.csv("/home/alphons/Uni/Forschungspraktikum/polar/scripts/u-content/uContentSmallTP.txt", header=FALSE)

bigTNdistribution = sort(bigUContentTN[, 1])
bigTPdistribution = sort(bigUContentTP[, 1])
smallTNdistribution = sort(smallUContentTN[, 1])
smallTPdistribution = sort(smallUContentTP[, 1])

par(mfrow=c(1,4))
hist(smallTNdistribution)
hist(smallTPdistribution)
hist(bigTNdistribution)
hist(bigTPdistribution)

#yTN = 1:length(TNdistribution)
#yTP = 1:length(TPdistribution)

#par(mfrow=c(1,2))
#plot(yTN, TNdistribution)
#plot(yTP, TPdistribution)