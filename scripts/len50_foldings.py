import sys
import tempfile
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


iput = open("../perf_testing/positiveSet.fa", "r")
seqLen = 50
offset = 260
oput = open("foldings_l" + str(seqLen) + "_o" + str(offset) + ".faf", "a")

def writeOutput(ipt, opt, length, off):
	faRecords = SeqIO.parse(ipt, "fasta")
	for record in faRecords:
		with tempfile.NamedTemporaryFile() as fp:
			#print(">" + record.id)
				
			fp.write(">" + record.id + "\n")
			#print(str(record.seq[off:off+length]))
			if record.id[-1] == "-":
				fp.write(str(record.seq[off+1:off+length+1].complement()))
			else:
				fp.write(str(record.seq[off:off+length]))
			fp.flush()
			cout = subprocess.check_output(["./RNAfold", "--noPS", "-i", fp.name])
			opt.write(cout)

if __name__ == "__main__":
	writeOutput(iput, oput, seqLen, offset)
