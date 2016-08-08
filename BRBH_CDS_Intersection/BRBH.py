import csv
import sys

"""
determine the best hit ortholgs after blast, please provide two four-colcum blast output, such as
ENSTRUT00000007201,ENST00000566487,194,5e-81
ENSTRUT00000007201,ENST00000422143,73.9,5e-81
ENSTRUT00000006789,ENST00000444775,54.7,5e-81
ENSTRUT00000026330,ENST00000390348,31.8,5e-81
ENSTRUT00000000001,ENST00000297142,29.5,5e-81
ENSTRUT00000000001,ENST00000297142,149,2e-41
ENSTRUT00000000001,ENST00000297142,35.0,2e-41
ENSTRUT00000000001,ENST00000297142,29.0,2e-41
ENSTRUT00000000001,ENST00000297142,97.3,2e-23
ENSTRUT00000000001,ENST00000297142,31.8,2e-23

and 
ENST00000390340,ENSTRUT00000004497,60.6,2e-09
ENST00000444775,ENSTRUT00000006789,66.9,3e-11
ENST00000390348,ENSTRUT00000026330,39.6,8e-07
ENST00000390348,ENSTRUT00000026330,29.5,8e-07
ENST00000566487,ENSTRUT00000007201,36.3,5e-07
ENST00000566487,ENSTRUT00000007201,25.8,5e-07
ENST00000566487,ENSTRUT00000007201,22.1,5e-07
ENST00000422143,ENSTRUT00000007201,36.3,5e-07
ENST00000422143,ENSTRUT00000007201,25.8,5e-07
ENST00000422143,ENSTRUT00000007201,22.1,5e-07

you would better firstly sorted your blast result by Evalue (12c) and score (11c)

"""


def blast_out_dict(blast_out):
	ortholog_dict = {}
	ortholog_listoflists = []
	for (query,subject,score,Evalue) in csv.reader(open(blast_out)):
		ortholog_list = []
		query = shorten_id(query)
		subject = shorten_id(subject)
		Evalue = float(Evalue)
		score = float(score)
		ortholog_list = [query,subject,score,Evalue]
		ortholog_listoflists.append(ortholog_list)
	# ortholog_listoflists_sort = sorted(sorted(ortholog_listoflists, key = lambda x: (x[0], x[2]), reverse = True), key = lambda x: (x[0], x[3]))
	ortholog_listoflists_sort = sorted(ortholog_listoflists, key = lambda x: (x[0], x[3], -x[2]))
	for lst in ortholog_listoflists_sort:
		query1,subject1,score1,Evalue2 = lst
		if query1 not in ortholog_dict and Evalue < 1e-3:
			ortholog_dict[query1] = subject1
	return ortholog_dict
			
		
		
"""
cut 'gi|XXXXX|' in NCBI format.
NCBI does this automatically for subject sequences.
"""		
def shorten_id(seqid):
    if seqid.startswith('gi|'):
        seqid = seqid.split('|')
        seqid = seqid[2:]
        seqid = "|".join(seqid)
    return seqid

	
"""find reciprocal best-hit"""	
# rbh = open("reciprocal_best_hit.txt",'w')
rbh = csv.writer(sys.stdout)
dict1 = blast_out_dict(sys.argv[1])
dict2 = blast_out_dict(sys.argv[2])

for id1 in dict1:
	# for id2 in dict2:
		# if dict1[id1] == id2 and dict2[id2] == id1:
			# rbh.writerow([id2, id1])

	match1 = dict1[id1]
	match2 = dict2.get(match1)
	if id1 == match2:
		rbh.writerow([id1,match1])












