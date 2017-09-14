import os
from os.path import join

def search_genes():
	matchingname=[]
	output=open("matching_genelines.txt",'w')
	genenames=open("genefunction.txt").readlines()
	genenames1=[name.replace("\n","") for name in genenames]
	print genenames1

	with open("ratebyGO10amph.txt") as file:
		allgenes=file.readlines()
		for genename in genenames1:
			for line in allgenes:
				if genename in line:
					output.write(line)
#							matchingname.append(genename)
#				print matchingname
if __name__ == "__main__":
    search_genes()