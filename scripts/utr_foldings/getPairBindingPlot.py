import os
import sys

def readBasePairProb(fileIter):
	bppDict = {}
	for line in fileIter:
		lineSplit = line.split(" ")
		try:
			iPos = int(lineSplit[0])
			jPos = int(lineSplit[1])
			prob = float(lineSplit[2])
		except ValueError:
			break
		
		ubox = lineSplit[3]
		if "ubox" in ubox:
			print (str(iPos) + " " + str(jPos) + " " + str(prob))
			bppDict[(iPos, jPos)] = prob
		else:
			break
	return bppDict

def readPolyaPosPerId(fileObj):
	ppDict = {}
	lineCount = 0
	for line in fileObj:
		if lineCount == 0:
			seqId = line[1:-1]
			lineCount += 1
			continue
		elif lineCount == 1: 
			polyaPosList = line[12:-2].split(",")
			lineCount += 1
			continue
		elif lineCount == 2:
			ppDict[seqId] = polyaPosList
			lineCount = 0
			continue
	return ppDict

def createCsvOutput(fileObj, bbpDict, pos, length):
	
	


if __name__ == "__main__":
	polyaF = open("../knownPolya.txt", "r")
	polyAPosDict = readPolyaPosPerId(polyaF)
	outputCsv = open("basePairProbAtPAS.csv", "a")

	for psFile in os.listdir("."):
		if psFile.endswith(".ps"):
			psF = open(psFile, "r")
			
			basePairProbDict = {}
			for line in f:
				if "%start of base pair probability data" in line:
					next(f)
					basePairProbDict = readBasePairProb(f)

			seqId = psFile.split("|")[0]
		
		for pasPos in polyAPosDict[seqId]
			createCsvOutput(outputCsv, basePairProbDict, pasPos, 100)	

