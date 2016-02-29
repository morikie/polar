from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import os

def countAdenin(seq):
	windowSize = 9
	aContentAtPos = []
	for i in range(0, len(seq) - windowSize):
		aContentAtPos.append(seq.count("a", i, i + windowSize))
	return aContentAtPos

#writes the Uracil content per sequence and position in a CSV file (overwrites the ouput file)
def writeCSV(fasta, outCSV, bigA, smallA, onlyCanonical):
	canonical = ["aataaa", "attaaa"]
	with open(outCSV, "w") as out, open(bigA, "w") as bigOut, open(smallA, "w") as smallOut:
		with open(fasta, "r") as f:
			lineCounter = 0
			negativeStrand = False
			for line in f:
				if lineCounter == 0:
					if line[-2] == "-":
						negativeStrand = True
					lineCounter += 1	
				elif lineCounter == 1:
					if negativeStrand:
						revComp = Seq(line, generic_dna)
						revComp = revComp.reverse_complement()
						if onlyCanonical and revComp[250:256] not in canonical:
							lineCounter = 0
							continue
						aContentList = countAdenin(str(revComp)[:-1])
						negativeStrand = False
						
						smallOut.write(str(revComp[256:265].count("a")))
						smallOut.write("\n")
						bigOut.write(str(revComp[260:300].count("a")))
						bigOut.write("\n")
					else:
						if onlyCanonical and line[250:256] not in canonical:
							lineCounter = 0
							continue
						aContentList = countAdenin(line[1:])
						
						smallOut.write(str(line[256:265].count("a")))
						smallOut.write("\n")
						bigOut.write(str(line[260:300].count("a")))
						bigOut.write("\n")

					out.write(str(aContentList)[1:-1])
					out.write("\n")

					lineCounter = 0
		
if __name__ == "__main__":
	tnSet = "../bin/TNdataSet.fa"
	tpSet = "../bin/TPdataSet.fa"
	fnSet = "../bin/FNdataSet.fa"
	fpSet = "../bin/FPdataSet.fa"
	negativeSet = "../perf_testing/negativeSet.fa"
	positiveSet = "../perf_testing/positiveSet.fa"

	resultCsvNeg = "resultAcontentNegative.csv"
	resultCsvPos = "resultAcontentPositive.csv"
	resultCsvTN = "resultAcontentTN.csv"
	resultCsvTP = "resultAcontentTP.csv"
	resultCsvFN = "resultAcontentFN.csv"
	resultCsvFP = "resultAcontentFP.csv"
	
	aContentBigTN = "aContentBigTN.txt"
	aContentBigTP = "aContentBigTP.txt"
	aContentSmallTN = "aContentSmallTN.txt"
	aContentSmallTP = "aContentSmallTP.txt"
	
	writeCSV(negativeSet, resultCsvNeg, aContentBigTN, aContentSmallTN, False)
	writeCSV(positiveSet, resultCsvPos, aContentBigTP, aContentSmallTP, False)
	writeCSV(tnSet, resultCsvTN, aContentBigTN, aContentSmallTN, False)
	writeCSV(tpSet, resultCsvTP, aContentBigTN, aContentSmallTN, False)
	writeCSV(fnSet, resultCsvFN, aContentBigTP, aContentSmallTP, False)
	writeCSV(fpSet, resultCsvFP, aContentBigTP, aContentSmallTP, False)



