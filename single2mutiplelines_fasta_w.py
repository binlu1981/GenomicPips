import sys
import os
inputfile = sys.argv[1]
length = int(sys.argv[2])

outfile = open(inputfile.replace(os.path.splitext(inputfile)[1],'_multi_line.fasta'),'w') #open outfile for writing

with open(inputfile, 'r') as f:
	for line in f:
		if line.startswith(">"):
			outfile.write(line)
		else:
			sequence = line.strip()
			while len(sequence) > 0:
				outfile.write(sequence[:length]+'\n')
				sequence = sequence[length:]