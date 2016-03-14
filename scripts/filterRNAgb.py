#author: Moritz Kiekeben
#email: moritz.k@fu-berlin.de
#description: 	
#Script to filter transcripts with known poly(A) signal sites from a GenBank file

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys


def writeToFile( path, id , seq, polyApos ):
	path.write(">" + id + "\n")
	path.write("Poly(A)_pos:")
	for pos in polyApos:
		path.write(str(pos) + ",")
	path.write("\n")
	path.write(str(seq) + "\n")


input_handle = open("rna.gbk", "rU")
output_handle = open("knownPolya.txt", "a")

total = 0
knownPolyA = 0

for record in SeqIO.parse(input_handle, "genbank"):
	total += 1
	polyApos = []
	
	for feat in record.features:
		regulatoryList = feat.qualifiers.get("regulatory_class")
		if record.id[0:3] == "NM_" and feat.type == "CDS":
			if feat.location and feat.location.end + 1 < len(record.seq):
				startPos = feat.location.start
				endPos = feat.location.end
				print (">" + record.id + "|UTR")
				print ((record.seq[endPos:]).lower())
		if regulatoryList != None:
			for item in regulatoryList:
				if item == "polyA_signal_sequence":
					polyApos.append(int(feat.location.start))
					knownPolyA += 1
	if polyApos:
		writeToFile(output_handle, record.id, record.seq, polyApos)
	



print ("total: ", total)
print ("knownPolyA: ", knownPolyA)

		
			
