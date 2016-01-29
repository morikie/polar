from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def countUracil(seq):
	windowSize = 9
	uContentAtPos = []
	for i in range(0, len(seq) - windowSize):
		uContentAtPos.append(seq.count("t", i, i + windowSize))
	return uContentAtPos

tnSet = "../perf_testing/tnSet.fa"
tpSet = "../perf_testing/tpSet.fa"
resultCsvTN = "resultUcontentTN.csv"
resultCsvTP = "resultUcontentTP.csv"

#TODO need fix: positions are no longer corresponding if the sequence is reverse complemented

with open(resultCsvTN, "a") as fResultTN, open(resultCsvTP, "a") as fResultTP:
	with open(tnSet, "r") as fTn:
		lineCounter = 0
		negativeStrand = False
		for line in fTn:
			if lineCounter == 0:
				if line[len(line)-1] == "-":
					negativeStrand = True
				lineCounter += 1	
			elif lineCounter == 1:
				if negativeStrand:
					revComp = Seq(line, generic_dna)
					revComp = revComp.reverse_complement()
					uContentList = countUracil(str(revComp))
					negativeStrand = False
				else:
					uContentList = countUracil(line)
				fResultTN.write(str(uContentList)[1:-1])
				fResultTN.write("\n")
				lineCounter = 0
	with open(tpSet, "r") as fTp:
		lineCounter = 0
		negativeStrand = False
		for line in fTp:
			if lineCounter == 0:
				if line[len(line)-1] == "-":
					negativeStrand = True
				lineCounter += 1	
			elif lineCounter == 1:
				if negativeStrand:
					revComp = Seq(line, generic_dna)
					revComp = revComp.reverse_complement()
					uContentList = countUracil(str(revComp))
					negativeStrand = False
				else:
					uContentList = countUracil(line)
				fResultTP.write(str(uContentList)[1:-1])
				fResultTP.write("\n")
				lineCounter = 0

