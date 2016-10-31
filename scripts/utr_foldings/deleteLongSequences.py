import sys
import os

iput = open("rna_UTRs.fa", "r")
oput = open("rna_UTRs_perfTrimmed.fa", "w")

# "0" being the header/sequence identifier and "1" being the sequence
lineIdentifier = 0
for line in iput:
    if lineIdentifier == 0:
        header = line
        lineIdentifier += 1
        continue
    if lineIdentifier == 1:
        if len(line) >= 2000:
            lineIdentifier = 0
            continue
        else:
            lineIdentifier = 0
            oput.write(header + line)
