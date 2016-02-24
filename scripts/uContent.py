from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import os

def countUracil(seq):
	windowSize = 9
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
	tnSet = "../perf_testing/tnSet.fa"
	tpSet = "../perf_testing/tpSet.fa"
	resultCsvTN = "resultUcontentTN.csv"
	resultCsvTP = "resultUcontentTP.csv"
	uContentBigTN = "uContentBigTN.txt"
	uContentBigTP = "uContentBigTP.txt"
	uContentSmallTN = "uContentSmallTN.txt"
	uContentSmallTP = "uContentSmallTP.txt"
	
	writeCSV(tnSet, resultCsvTN, uContentBigTN, uContentSmallTN, False)
	writeCSV(tpSet, resultCsvTP, uContentBigTP, uContentSmallTP, False)
