from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def writeFilteredFasta(inFasta, outFasta):
	canonical = ["aataaa", "attaaa"]
	with open(inFasta, "r") as inFa, open(outFasta, "w") as outFa:
		headerLine = ""
		for line in inFa:
			if line[0] == ">":
				headerLine = line
			else:
				strand = headerLine[-2]
				if strand == "-":
					revComp = Seq(line, generic_dna)
					line = str(revComp.reverse_complement())[1:-1]
				else:
					line = line[1:-1]
				#print(strand + ": " + line)
				if line[249:255] in canonical:
					outFa.write(headerLine + line + "\n")

if __name__ == "__main__":
	tpSet = "../perf_testing/tpSet.fa"
	tnSet = "../perf_testing/tnSet.fa"
	canonicalTpSet = "../perf_testing/canonicalTpSet.fa"
	canonicalTnSet = "../perf_testing/canonicalTnSet.fa"

	writeFilteredFasta(tpSet, canonicalTpSet)
	writeFilteredFasta(tnSet, canonicalTnSet)
