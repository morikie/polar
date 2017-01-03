dataFolder = "/home/alphons/Uni/polar/scripts/"
tempPos = read.csv(paste0(dataFolder, "utr_foldings/basePairProbAtPAS.csv"), header=FALSE)

xPos = seq(0, 99)
yPos = rep(0, 100)

for (i in 1:(length(yPos) - 1)) {
  yPos[i] = mean(tempPos[, i])
}
par(mfrow=c(1,1), oma=c(0,0,4,0))
plot(xPos, 
     yPos, 
     main="Base pair probability data set (average)", 
     type="h", 
     ylim=c(0,0.5))

axis(side=1, at=seq(0, 100, by=5))
abline(v=75, col="red")
abline(v=82, col="red")

motif = seq(76, 81)
#motif = seq(1, 6)
par(mfrow=c(3,2), oma=c(0,0,4,0))
for (i in 1:length(motif)) {
  hist(tempPos[,motif[i]], 
       ylim=c(0, 1500), 
       breaks=seq(0, 1, by=0.002), 
       xlim=c(0.0, 0.2), 
       main=paste0("BPP distribution at positon ", motif[i]), 
       xlab="BPP")
}
mtext("BPP distributions at corresponding motif positions", outer=TRUE, cex=1.5)

par(mfrow=c(1,1), oma=c(0,0,4,0))
motif = seq(76, 81)
sumMotifBpp = {}
for (i in 1:length(tempPos[, 1])) {
  sumMotifBpp = append(sumMotifBpp, sum(tempPos[i, motif]))
}

hist(sumMotifBpp,
     ylim=c(0,2500),
     breaks=seq(0,8,by=0.1), 
     main="Histogram of sums of BPPs of 6 consecutive positions in each sequence", 
     xlab="sum of 6 consec. BPPs")
