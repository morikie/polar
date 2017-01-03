dataFolder = "/home/alphons/Uni/polar/"
truthValues = read.csv(paste0(dataFolder, "bin/sensitivityOut.csv"), header=FALSE)

hist(truthValues[,2])

