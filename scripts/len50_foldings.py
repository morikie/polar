import sys
import tempfile
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


iput = open("../bin/TPdataSet.fa", "r")
oput = open("plus20_len50_foldings.faf", "a")
seqLen = 50
offset = 276

def writeOutput(ipt, opt, length, off):
	faRecords = SeqIO.parse(ipt, "fasta")
	for record in faRecords:
		with tempfile.NamedTemporaryFile() as fp:
			#print(">" + record.id)
			fp.write(">" + record.id + "\n")
			#print(str(record.seq[off:off+length]))
			fp.write(str(record.seq[off:off+length]))
			fp.flush()
			cout = subprocess.check_output(["./RNAfold", "-i", fp.name])
			opt.write(cout)

if __name__ == "__main__":
	writeOutput(iput, oput, seqLen, offset)
