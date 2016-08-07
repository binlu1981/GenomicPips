import glob
import csv
import os

"""
this script was used to extract shared ortholgs from each fasta dataset in all subfolders, 
then merge them according reference id.
Users should prepare 1vs1_ortholgos.txt as csv format in each subfolders, like as:

aaa,ENST00000000233
ade,ENST00000000412
bt,ENST00000001008

and fasta set in each subfolders as , like as:
>aaa
atcgagtcaaa
>ade
cgtcggaaccc
>bt
gtggtcccccagt

For example, thare are three subfolders represent three different species, 
we plan to find the shared genes among them, and then merge as individual fasta file.
file distribution like this:

ac	gt	kli	sharefilesmergelist.py

./ac:
ac_ortholog.fas  ac1vs1_orthologs.txt

./gt:
gt_ortholog.fas  gt1vs1_orthologs.txt

./kli:
kli_ortholog.fas  kli1vs1_orthologs.txt

just to do
python sharefilesmergelist.py
"""



folderpathdict = {}
# folderiddict = {}
folderseqdict = {}
folderrefiddictdict = {}
# idseqdict = {}
# seqsetlist = []
# idreflist = []
refiddictlist = []
cwd = os.getcwd()
folderlist = next(os.walk(cwd))[1]
for folder in folderlist:
	folderpath = os.path.join(cwd,folder)
	filelist = os.listdir(folderpath)
	folderpathdict[folder] = folderpath
	for file in filelist:		
		if file.endswith("orthologs.txt"):
			# idreflist.append(file)
			# folderiddict[folder] = file
			refiddict = {}
			for (id,ref) in csv.reader(open(os.path.join(cwd,folder,file))):
				refiddict[ref] = id
			refiddictlist.append(refiddict)
			folderrefiddictdict[folder] = refiddict
		if file.endswith("ortholog.fas"):
			# seqsetlist.append(file)
			folderseqdict[folder] = file
		# idseqdict = dict(zip(idreflist,seqsetlist))

reflist = [list(d.keys()) for d in refiddictlist]
shareref = set.intersection(*map(set,reflist))
# shareref = reduce(set.intersection, map(set, [list(d.values()) for d in refiddictlist]))
# shareref = reduce(lambda x, y: x&y, (set(d.values()) for d in refiddictlist))

for sharerefid in shareref:
	shareseqcombine = open(sharerefid+".fas","w")
	for ws in folderseqdict:
		seq1 = open(os.path.join(folderpathdict[ws],folderseqdict[ws])).read()
		seq2 = seq1.split(">")
		seq3 = seq2[1:]
		seq4 = [seq.partition("\n") for seq in seq3]
		seq5 = {seq[0]:seq[2].replace("\n","") for seq in seq4}
		# seq5 = [[seq[0],seq[2].replace("\n","")] for seq in seq4]
		# seq6 = dict((seq[0],seq[1]) for seq in seq5)	

		if sharerefid in folderrefiddictdict[ws].keys():
			seqname = folderrefiddictdict[ws][sharerefid]
			shareseq = seq5[seqname]
		shareseqcombine.write("%s\n%s\n" % (">"+seqname,shareseq))


	
	

	
# print folderrefiddictdict
# print reflist
# print shareref		
# print refiddictlist
# print folderpathdict
# print folderseqdict


	
	
	