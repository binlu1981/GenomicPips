import re
def flat(tree):
    res = []
    for i in tree:
        if isinstance(i, list):
            res.extend(flat(i))
        else:
            res.append(i)
    return res
	

def read_davidcluster():
	print "input your file name"
	t2t=raw_input(">")
	print "input your interesting cluster number"
	number=raw_input(">")
	output=open("david_cluster"+number+"_genes.txt",'w')
	with open(t2t) as file:
		txt=file.read()
#		print txt
		txt1=txt.split("Annotation ")
#		print txt1
		header=("Cluster "+number+"\t")
		geneid=("ENST.*")
		for txt2 in txt1:
			if header in txt2:
#				print txt2
				line0=re.findall("ENST.*",txt2)
#				print line0
				line1=[line.partition("\t") for line in line0]
#				print line1
				line2=[line[0] for line in line1]
#				print line2
				line3=list(line.split(",") for line in line2)
#				print line3
				line4=list(set(flat(line3)))
				print line4
				for gene in line4:
					geneid=''.join(gene.split())
					print geneid
					output.write(geneid+"\n")
#				line4=str(line3).replace("[",'').replace("]",'')
#				print line4
#				print line4.split(',')
#				return [eval(x) for x in line4.split(',') if x.strip()]
#				gene=list(set(line2))
#				print gene
#				txt3=txt2.split("\n")
#				print txt3
#				if geneid in txt3:
#					print geneid
				
				
if __name__ == "__main__":
    read_davidcluster()
