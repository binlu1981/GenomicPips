import re

tn=raw_input("enrichment file name? ")
d=raw_input("interesting Term name? ")

id=open(d)
names=id.readlines()
names1=[n.strip() for n in names]
print names1

interesting=open("interestingenrichment_DAVID.txt",'w')

with open(tn) as en:
	lst= en.read().split("\n\n")
	lst1=filter(None, lst)
	patternscore=re.compile("Enrichment\sScore\:\s+(\d+\.\d+)\n")
	patternineachcluster=re.compile("FDR\n(.+)",re.DOTALL)
	for lst2 in lst1:
		score=re.findall(patternscore,lst2)
		cluster=re.findall(patternineachcluster,lst2)
#		print cluster
		eachline=cluster[0].split("\n")
		line1=[l+'\t'+''.join(score) for l in eachline]
		for l in line1:
			line1a=''.join(l)
			for name in names1:
				if name in line1a:
					print line1a
					interesting.write("%s\n" % (line1a))
