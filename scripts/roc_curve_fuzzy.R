library(ROCR)
temp = read.csv("/home/morikie/Documents/Forschungspraktikum/C++/polar/bin/fuzzyOut.txt", header=TRUE)
senVec = temp[,1]
fprVec = 1 - temp[,2]
predVec = append(senVec, fprVec)
labelVec = append(rep(1, length(senVec)), rep(0, length(fprVec)))
pred = prediction(predVec, labelVec)
roc = performance(pred, "tpr", "fpr")
auc = performance(pred, measure="auc")
plot(roc)
auc@y.values