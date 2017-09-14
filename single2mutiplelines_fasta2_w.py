import sys
import os
inputfile = sys.argv[1]
length = int(sys.argv[2])

print '\nConverting ' + inputfile + ' to ' + str(length) + ' chars per line...'

with open(inputfile, 'r') as f:

	h = []                                             #headers list
	s = []                                             #sequence list
	for line in f:
		if line.startswith(">"):                      #grab headers
			h.append(line.strip().split(">")[1])
		else:
			s.append(line.strip())                   #grab sequences

myseqs = dict(zip(h,s))                                 #map to dictionary, header:sequence

outfile = open(inputfile.replace(os.path.splitext(inputfile)[1],'_multi_line.fasta'),'w') #open outfile for writing

def print_seqs(header, sequence):
		outfile.write('>'+header+'\n')
		while len(sequence) > 0:
			outfile.write(sequence[:length]+'\n')
			sequence = sequence[length:]

for i in myseqs:
		print_seqs(i, myseqs[i])
print 'Converted ' + str(len(myseqs)) + ' sequences.'
print 'Done.\n'

outfile.close()