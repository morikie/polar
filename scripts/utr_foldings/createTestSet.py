f = open("rna_UTRs.fa", "r")
fout = open("rna_small_UTRs.fa", "w")

i = 0
prevLine = ""
matches = 0
for line in f:
	if i == 0:
		prevLine = line
		i += 1
		continue
	if i == 1:
		if line.find("aataaa", len(line) - 50) == -1:
			i = 0
			continue
		if matches > 20:
			break
		fout.write(prevLine)
		fout.write(line)
		matches += 1
		i = 0	
f.close()
fout.close()
	
