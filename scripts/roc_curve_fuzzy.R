library(ROCR)
temp = read.csv("/home/alphons/Uni/Forschungspraktikum/polar/bin/perfOut_newFeature.txt", header=TRUE)
senVec = temp[,2]
fprVec = 1 - temp[,3]
predVec = append(senVec, fprVec)
labelVec = append(rep(1, length(senVec)), rep(0, length(fprVec)))
pred = prediction(predVec, labelVec)
roc = performance(pred, "sens", "fpr")
auc = performance(pred, measure="auc")
par(mfrow=c(1,1))
plot(roc)
abline(0,1)
auc@y.values