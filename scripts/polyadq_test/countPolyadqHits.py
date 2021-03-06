
def relativeHits(f, hitPos):
	hits = 0
	numSeq = 0
	with open(f, "r") as f:
		foundPositive = False
		for num, line in enumerate(f, 1):
			if line[0:3] == "chr":
				numSeq += 1
			if "<td align=center>POS" in line:
				foundPositive = True
				continue
			if "<td align=center>" + str(hitPos) in line and foundPositive:
				hits += 1
			foundPositive = False
	print(hits)
	print(numSeq)
	print("")
	return float(hits)/numSeq

if __name__ == "__main__":
	polyadqTnHtml = "polyadq4-TN.cgi.html"
	polyadqTpHtml = "polyadq4-TP.cgi.html"
	
	print("checking positive set")
	sensitivity = relativeHits(polyadqTpHtml, 250)
	print("checking negative set")
	specificity = 1 - relativeHits(polyadqTnHtml, 251)
	
	print("sensitivity: " + str(sensitivity))
	print("specificity: " + str(specificity))

