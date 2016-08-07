import os
import csv
import glob

"""extract cds region accord position information in tblastx output"""

seqset = ''.join(glob.glob("*ortholog.fas"))
blastinfo = ''.join(glob.glob("*.tab_sortU"))

cdsset = open(seqset.replace(".fas","_cds.fas"),'w')

for (query_id, subject_id, identity, alignment_length, mismatches, gap_opens, q_start, q_end, s_start, s_end, evalue, bit_score) in csv.reader(open(blastinfo),delimiter='\t'):
	siminfo = [query_id, subject_id, q_start, q_end]
	# print siminfo
	seq1 = open(seqset).read()
	seq2 = seq1.split(">")
	seq3 = seq2[1:]
	seq4 = [seq.partition("\n") for seq in seq3]
	seq5 = {seq[0]:seq[2].replace("\n","") for seq in seq4}
	# seq5 = [[seq[0],seq[2].replace("\n","")] for seq in seq4]
	# seq6 = dict((seq[0],seq[1]) for seq in seq5)
	for gene in seq5.keys():
		if gene == query_id:
			cdsind = open(gene+".fas",'w')
			cdna = seq5[gene]
			if q_end >= q_start:
				cds = cdna[int(q_start)-1:int(q_end)]
			else:
				cds_reverse_complement = cdna[int(q_end)-1:int(q_start)]
				complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N': 'N','X': 'X','U': 'A','-': '-','a': 't', 'c': 'g', 'g': 'c', 't': 'a','n': 'n','x': 'x','u': 'a'}
				cds = "".join(complement.get(base, base) for base in reversed(cds_reverse_complement))

			print gene, cds
			cdsset.write("%s\n%s\n" % (">"+gene,cds))
			cdsind.write("%s\n%s\n" % (">"+gene,cds))


	






