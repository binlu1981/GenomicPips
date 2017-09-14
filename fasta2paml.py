import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-e',"--extension",type=str,default='fas',help='enter file extension')
args = parser.parse_args()

print("USAGE: python fasta2paml.py -e <input file extension>") 

for file in glob.glob("*"+args.extension):
	outfile = file.replace(args.extension,"phylip")
	with open(file) as f:
		records = f.read()
		seq_list = records.split(">")[1:]
		seq1 = [seq.partition("\n") for seq in seq_list]
		seq2 = {'_'.join(seq[0].strip().split()):seq[2].replace('\n','').strip() for seq in seq1}
	with open(outfile,'w') as op:
#		op.write("\t%s  %s\n" % (len(seq2),len(seq2.values()[0])-1))
		op.write("\t{}  {}\n".format(len(seq2),len(list(seq2.values())[0])-1))
		for x in seq2:
#			op.write("%s  %s\n" % (x.strip(),seq2[x].strip()))
#			op.write("{:20}{}\n".format(x.strip(),seq2[x].strip()))
			op.write("{}  {}\n".format(x.strip(),seq2[x].strip()))
	
	