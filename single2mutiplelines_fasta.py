import sys

inputfile = sys.argv[1]
length = int(sys.argv[2])

outfile = open(inputfile.split(".fasta")[0] + '_multi-line.fasta', 'w') #open outfile for writing

with open(inputfile, 'r') as f:
	for line in f:
		if line.startswith(">"):
			print >> outfile, line.strip()
		else:
			sequence = line.strip()
			while len(sequence) > 0:
				print >>outfile, sequence[:length]
				sequence = sequence[length:]