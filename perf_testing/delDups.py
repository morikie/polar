import os
import sys

fIn = open("bppDict.txt", "r")
fOut = open("utrBppPerTranscriptNoDupes.txt", "w")

bppDict = {}

lineCounter = 0
for line in fIn:
    if lineCounter == 0:
        seqId = line[1:]
        lineCounter = 1
        continue
    if lineCounter == 1:
        if seqId in bppDict:
            lineCounter = 0
            seqId = ""
            continue
        else:
            bppDict[seqId] = line
            lineCounter = 0
            continue

for key, value in bppDict.items():
    fOut.write(">" + key + value)

