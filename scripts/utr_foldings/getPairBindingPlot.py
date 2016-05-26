import os
import sys
from collections import OrderedDict

class ijProbability:
	def __init__(self, i):
		self.i = i
		self.jList = []
		self.probList = []
	
	def add_j_prob(self, j, prob):
		self.jList.append(j)
		self.probList.append(prob)

#extracting the base pair probabilities from the postscript files that Vienna created
def readBasePairProb(fileIter):
	bppDict = OrderedDict()
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
			#print (str(iPos) + " " + str(jPos) + " " + str(prob))
			if iPos not in bppDict:
				ijProb = ijProbability(iPos)
				ijProb.add_j_prob(jPos, prob)
				bppDict[iPos] = ijProb
			else:
				bppDict[iPos].add_j_prob(jPos, prob)
		else:
			break
	return bppDict

#extracting PAS positions for each UTR sequence
def readPolyaPosPerId(fileObj):
	ppDict = {}
	lineCount = 0
	for line in fileObj:
		if lineCount == 0:
			seqId = line[1:-1]
			lineCount += 1
			continue
		elif lineCount == 1: 
			temp = line[12:-2].split(",")
			polyaPosList = []
			for pos in temp:
				try:
					polyaPosList.append(int(pos))
				except ValueError:
					print("ValueError: polyaPosList")
			lineCount += 1
			continue
		elif lineCount == 2:
			ppDict[seqId] = polyaPosList
			lineCount = 0
			continue
	return ppDict

def createCsvOutput(fileObj, bppDict, pos, length):
	bppSize = len(bppDict)
	outputList = []
	start = pos - length + 6
	print(str(pos) + ", " + str(start))
	for i in range(start, start + 100):
		if i in bppDict:
			outputList.append(max(bppDict[i].probList))
		else:
			outputList.append(0.0)
	
	fileObj.write(str(outputList)[1:-1] + "\n")

if __name__ == "__main__":
	polyaF = open("../knownPolya.txt", "r")
	polyAPosDict = readPolyaPosPerId(polyaF)
	outputCsv = open("basePairProbAtPAS.csv", "a")

	for psFile in os.listdir("."):
		if psFile.endswith(".ps"):
			psF = open(psFile, "r")
			
			basePairProbDict = []
			for line in psF:
				if "%start of base pair probability data" in line:
					next(psF)
					basePairProbDict = readBasePairProb(psF)
			headerSplit = psFile.split("|")
			seqId = headerSplit 
			offset = psFile.split("|")[2]
			offset = int(offset[7:-6])
		if seqId in polyAPosDict:	
			for pasPos in polyAPosDict[seqId]:
				sys.stdout.write(seqId + ", ")
				createCsvOutput(outputCsv, basePairProbDict, pasPos - offset, 100)	

