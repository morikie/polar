from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import os

def countUracil(seq):
	windowSize = 10
	uContentAtPos = []
	for i in range(0, len(seq) - windowSize):
		uContentAtPos.append(seq.count("t", i, i + windowSize))
	return uContentAtPos

#writes the Uracil content per sequence and position in a CSV file (overwrites the ouput file)
def writeCSV(fasta, outCSV, bigU, smallU, onlyCanonical):
	canonical = ["aataaa", "attaaa"]
	with open(outCSV, "w") as out, open(bigU, "w") as bigOut, open(smallU, "w") as smallOut:
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
						uContentList = countUracil(str(revComp)[:-1])
						negativeStrand = False
						
						smallOut.write(str(revComp[256:265].count("t")))
						smallOut.write("\n")
						bigOut.write(str(revComp[260:300].count("t")))
						bigOut.write("\n")
					else:
						if onlyCanonical and line[250:256] not in canonical:
							lineCounter = 0
							continue
						uContentList = countUracil(line[1:])
						
						smallOut.write(str(line[256:265].count("t")))
						smallOut.write("\n")
						bigOut.write(str(line[260:300].count("t")))
						bigOut.write("\n")

					out.write(str(uContentList)[1:-1])
					out.write("\n")

					lineCounter = 0
		
if __name__ == "__main__":
	tnSet = "../bin/TNdataSet.fa"
	tpSet = "../bin/TPdataSet.fa"
	fnSet = "../bin/FNdataSet.fa"
	fpSet = "../bin/FPdataSet.fa"
	negativeSet = "../perf_testing/negativeSet.fa"
	positiveSet = "../perf_testing/positiveSet.fa"

	resultCsvNeg = "u-content/resultUcontentNegative.csv"
	resultCsvPos = "u-content/resultUcontentPositive.csv"
	resultCsvTN = "u-content/resultUcontentTN.csv"
	resultCsvTP = "u-content/resultUcontentTP.csv"
	resultCsvFN = "u-content/resultUcontentFN.csv"
	resultCsvFP = "u-content/resultUcontentFP.csv"
	
	uContentBigTN = "u-content/uContentBigTN.txt"
	uContentBigTP = "u-content/uContentBigTP.txt"
	uContentSmallTN = "u-content/uContentSmallTN.txt"
	uContentSmallTP = "u-content/uContentSmallTP.txt"
	
	writeCSV(negativeSet, resultCsvNeg, uContentBigTN, uContentSmallTN, False)
	writeCSV(positiveSet, resultCsvPos, uContentBigTP, uContentSmallTP, False)
	writeCSV(tnSet, resultCsvTN, uContentBigTN, uContentSmallTN, False)
	writeCSV(tpSet, resultCsvTP, uContentBigTN, uContentSmallTN, False)
	writeCSV(fnSet, resultCsvFN, uContentBigTP, uContentSmallTP, False)
	writeCSV(fpSet, resultCsvFP, uContentBigTP, uContentSmallTP, False)



