#Example 2-16.
"""
from random import randint
def base_random_replace(seq):
	position=randint(0,len(seq)-1)
	seq1=seq[0:position]+"ATCG".replace(seq[position],"")[randint(0,2)]+seq[position+1:]
	return seq1.upper()
print base_random_replace("atcctggt")

#3-2 3-4
print set("atcg")
DNAbase={"a","t","c","g"}
RNAbase=set("aucg")
print DNAbase
print DNAbase|RNAbase
print DNAbase&RNAbase
print DNAbase^RNAbase
print DNAbase-RNAbase
print DNAbase.update(DNAbase|RNAbase) #DNAbase&|RNAbase
print DNAbase.intersection_update(DNAbase&RNAbase) #DNAbase&=RNAbase
print DNAbase.difference_update(DNAbase-RNAbase) #DNAbase-=RNAbase
print DNAbase.symmetric_difference_update(DNAbase^RNAbase)#DNAbase^=RNAbase
#Example 3-2
def restriction_cut(seq,recognition_seq):
	position=seq.find(recognition_seq)
	return seq[0:position],seq[position:]
print restriction_cut("atcctggt","ctg")
#table 3-12
lst=["ac","ag"]
lst.append("gc")
print lst #["ac","ag","gc"]
lst.extend("gc") #["ac","ag","gc","g","c"]
print lst
lst.insert(1,"gc")
print lst #["ac","gc","ag","gc","g","c"]
lst.remove("gc")
print lst #["ac","ag","gc","g","c"] first "gc"
lst.pop(3)
print lst #["ac","ag","gc","c"]
lst.reverse()
print lst
lst.sort()
print lst

print "a\nc".splitlines(1)
print ">".join("ab")
print ">a\n>b".split(">") #['','a\n','b']
print ">a\n>b".rsplit(">",1) #['>a\n','b']
print ">a\n>b".partition(">") #('','>','a\n>b')
print ">a\n>b".rpartition(">") #('>a\n','>','b') reverse partition

#example3-5
RNA_codon_table = {
# Second Base
# U C A G
# U
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys', # UxU
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys', # UxC
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---', # UxA
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Urp', # UxG
# C
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg', # CxU
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg', # CxC
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg', # CxA
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg', # CxG
# A
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser', # AxU
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser', # AxC
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg', # AxA
'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg', # AxG
# G
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly', # GxU
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly', # GxC
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly', # GxA
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly' # GxG
}

print RNA_codon_table
from pprint import pprint
pprint(RNA_codon_table)

print open("fasta.fas",'r').readline()
print open("fasta.fas",'r').readlines()
print open("fasta.fas",'r').read()
#example 3-6
def read_fasta_strings(fastafile):
	with open(fastafile) as file:
		return file.read().split(">")[1:]
print read_fasta_strings("fasta.fas")
#['gi|ddfd566\naaaaaaaaaaaaaaaccccccccccccttttttttttttttggggggggggg\n', 'abs456\ntttttttttccccccccccccgggggggggggg\n', '1234\nttttttttccccccccccccccccccgggggggggggggggggaaaaaaaaaaa']
def read_fast_entries(fastafile):
	return [seq.partition('\n') for seq in read_fasta_strings(fastafile)]
print read_fast_entries("fasta.fas")
#[('gi|ddfd566', '\n', 'aaaaaaaaaaaaaaaccccccccccccttttttttttttttggggggggggg\n'), ('abs456', '\n', 'tttttttttccccccccccccgggggggggggg\n'), ('1234', '\n', 'ttttttttccccccccccccccccccgggggggggggggggggaaaaaaaaaaa')]
def read_fasta_sequences(fastafile):
	return [[seq[0],seq[2].replace('\n','')] for seq in read_fast_entries(fastafile)]
print read_fasta_sequences("fasta.fas")
def read_fasta_sequences1(fastafile):
	return [(info,seq.replace('\n','')) for info,ignore,seq in read_fast_entries(fastafile)]
print read_fasta_sequences1("fasta.fas")
#[('gi|ddfd566', 'aaaaaaaaaaaaaaaccccccccccccttttttttttttttggggggggggg'), ('abs456', 'tttttttttccccccccccccgggggggggggg'), ('1234', 'ttttttttccccccccccccccccccgggggggggggggggggaaaaaaaaaaa')]
def read_fasta_and_info(fastafile):
	return [[seq[0].split('|'),seq[1]] for seq in read_fasta_sequences(fastafile)]
print read_fasta_and_info("fasta.fas")
#[[['gi', 'ddfd566'], 'aaaaaaaaaaaaaaaccccccccccccttttttttttttttggggggggggg'], [['abs456'], 'tttttttttccccccccccccgggggggggggg'], [['1234'], 'ttttttttccccccccccccccccccgggggggggggggggggaaaaaaaaaaa']]

#example 3-15
def read_seq(fasta_file):
	seq=open(fasta_file).read()
	print "seq is \n%s" % seq
	seq1=seq.split(">")[1:]
	print "seq1 is \n%s" % seq1
	seq2=[seqs.partition("\n") for seqs in seq1]
	print "seq2 is \n%s" % seq2
	seq3=[(seqs[0],seqs[2].replace('\n','')) for seqs in seq2]
	print "seq3 is \n%s" % seq3
	seq4=[(seqs[0].split('|'),seqs[1]) for seqs in seq3]
	print "seq4 is \n%s" % seq4
read_seq("fasta.fas")
#seq is
#>gi|ddfd566
#aaaaaaaaaaaaaaaccccccccccccttttttttttttttggggggggggg
#>abs456
#tttttttttccccccccccccgggggggggggg
#>1234
#ttttttttccccccccccccccccccgggggggggggggggggaaaaaaaaaaa
#seq1 is
#['gi|ddfd566\naaaaaaaaaaaaaaaccccccccccccttttttttttttttggggggggggg\n', 'abs456\ntttttttttccccccccccccgggggggggggg\n', '1234\nttttttttccccccccccccccccccgggggggggggggggggaaaaaaaaaaa']
#seq2 is
#[('gi|ddfd566', '\n', 'aaaaaaaaaaaaaaaccccccccccccttttttttttttttggggggggggg\n'), ('abs456', '\n', 'tttttttttccccccccccccgggggggggggg\n'), ('1234', '\n', 'ttttttttccccccccccccccccccgggggggggggggggggaaaaaaaaaaa')]
#seq3 is
#[('gi|ddfd566', 'aaaaaaaaaaaaaaaccccccccccccttttttttttttttggggggggggg'), ('abs456', 'tttttttttccccccccccccgggggggggggg'), ('1234', 'ttttttttccccccccccccccccccgggggggggggggggggaaaaaaaaaaa')]
#seq4 is
#[(['gi', 'ddfd566'], 'aaaaaaaaaaaaaaaccccccccccccttttttttttttttggggggggggg'), (['abs456'], 'tttttttttccccccccccccgggggggggggg'), (['1234'], 'ttttttttccccccccccccccccccgggggggggggggggggaaaaaaaaaaa')]

#example 3-7
def validate_base(seq,RNA=False):
	valid="AUCG" if RNA else "ATCG"
	return all([(base in valid) for base in seq.upper()])
print validate_base("atcccgggggttt")
#example 3-8
from random import randint
def random_base(RNA=False):
	return ("AUCG" if RNA else "ATCG")[randint(0,3)]
print random_base()
def random_codon(RNA=False):
	return random_base(RNA)+random_base(RNA)+random_base(RNA)
print random_codon(RNA=False)
print [random_codon() for n in range(randint(3,10))]
def random_codons(minlen=3,maxlen=10,RNA=False):
	return [random_codon(RNA) for n in range(randint(minlen,maxlen))]
print random_codons()
#example 3-9
def translate_RNA_codon(codon):
	return RNA_codon_table[codon]
#example 3-9
def random_codon_translation(minlen=3,maxlen=10):
	return [RNA_codon_table[codon] for codon in random_codons(minlen,maxlen,True)]
print random_codon_translation()

#example 3-17
def aa_generate(seq):
	for base in seq:
		if base in set("AUCGaugc"):
			return [RNA_codon_table[seq.upper()[n:n+3]] for n in range(0,len(seq),3)]
		else:
			seq1=seq.upper().replace("T","U")
			return [RNA_codon_table[seq1[n:n+3]] for n in range(0,len(seq),3)]
aagen=aa_generate("tcACCGCACCAACAGCGC")
print aagen



#example 3-19
def get_seq_description(fastafile):
	seq=open(fastafile)
	header=[]
	for line in seq:		
		if line[0]==">":
			description=line[1:].split("|")
			header.append(description)
	return header
print get_seq_description("fasta.fas")
def get_seq_description1(fastafile):
	seq=open(fastafile)
	return [line[1:].split("|") for line in seq if line[0]==">"]
print get_seq_description1("fasta.fas")

def first_common(list1, list2):
	return [item for item in list1 if item in list2]
list1=[1,2,3,67]
list2=[1,3,6,7]
print first_common(list1, list2)

def get_dna_code(fastfile):
	with open(fastfile) as file:
		for line in file:
			if len(line.split('|'))<1:
				return []
		return {line[2] for line in file if line[0] == '>'}
print get_dna_code("fasta.fas")
			
#example 3-23
def generate_triples(chars='TCAG'):
	charss=set(chars)
	return [b1+b2+b3 for b1 in charss for b2 in charss for b3 in chars]
print generate_triples()
print generate_triples("01")

print max(range(3,7),key=abs)

def fn(x,y):
	return x*x+y*y
fn1=lambda x,y: x*x+y*y
print fn(5,6)
print fn1(5,6)
"""
standardWithAmbig = {
	'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAR':'K', 'AAT':'N', 'AAY':'N', 'ACA':'T', 'ACB':'T', 
	'ACC':'T', 'ACD':'T', 'ACG':'T', 'ACH':'T', 'ACK':'T', 'ACM':'T', 'ACN':'T', 'ACR':'T', 
	'ACS':'T', 'ACT':'T', 'ACV':'T', 'ACW':'T', 'ACY':'T', 'AGA':'R', 'AGC':'S', 'AGG':'R', 
	'AGR':'R', 'AGT':'S', 'AGY':'S', 'ATA':'I', 'ATC':'I', 'ATG':'M', 'ATH':'I', 'ATM':'I', 
	'ATT':'I', 'ATW':'I', 'ATY':'I', 'CAA':'Q', 'CAC':'H', 'CAG':'Q', 'CAR':'Q', 'CAT':'H', 
	'CAY':'H', 'CCA':'P', 'CCB':'P', 'CCC':'P', 'CCD':'P', 'CCG':'P', 'CCH':'P', 'CCK':'P', 
	'CCM':'P', 'CCN':'P', 'CCR':'P', 'CCS':'P', 'CCT':'P', 'CCV':'P', 'CCW':'P', 'CCY':'P', 
	'CGA':'R', 'CGB':'R', 'CGC':'R', 'CGD':'R', 'CGG':'R', 'CGH':'R', 'CGK':'R', 'CGM':'R', 
	'CGN':'R', 'CGR':'R', 'CGS':'R', 'CGT':'R', 'CGV':'R', 'CGW':'R', 'CGY':'R', 'CTA':'L', 
	'CTB':'L', 'CTC':'L', 'CTD':'L', 'CTG':'L', 'CTH':'L', 'CTK':'L', 'CTM':'L', 'CTN':'L', 
	'CTR':'L', 'CTS':'L', 'CTT':'L', 'CTV':'L', 'CTW':'L', 'CTY':'L', 'GAA':'E', 'GAC':'D', 
	'GAG':'E', 'GAR':'E', 'GAT':'D', 'GAY':'D', 'GCA':'A', 'GCB':'A', 'GCC':'A', 'GCD':'A', 
	'GCG':'A', 'GCH':'A', 'GCK':'A', 'GCM':'A', 'GCN':'A', 'GCR':'A', 'GCS':'A', 'GCT':'A', 
	'GCV':'A', 'GCW':'A', 'GCY':'A', 'GGA':'G', 'GGB':'G', 'GGC':'G', 'GGD':'G', 'GGG':'G', 
	'GGH':'G', 'GGK':'G', 'GGM':'G', 'GGN':'G', 'GGR':'G', 'GGS':'G', 'GGT':'G', 'GGV':'G', 
	'GGW':'G', 'GGY':'G', 'GTA':'V', 'GTB':'V', 'GTC':'V', 'GTD':'V', 'GTG':'V', 'GTH':'V', 
	'GTK':'V', 'GTM':'V', 'GTN':'V', 'GTR':'V', 'GTS':'V', 'GTT':'V', 'GTV':'V', 'GTW':'V', 
	'GTY':'V', 'MGA':'R', 'MGG':'R', 'MGR':'R', 'TAA':'*', 'TAC':'Y', 'TAG':'*', 'TAR':'*', 
	'TAT':'Y', 'TAY':'Y', 'TCA':'S', 'TCB':'S', 'TCC':'S', 'TCD':'S', 'TCG':'S', 'TCH':'S', 
	'TCK':'S', 'TCM':'S', 'TCN':'S', 'TCR':'S', 'TCS':'S', 'TCT':'S', 'TCV':'S', 'TCW':'S', 
	'TCY':'S', 'TGA':'*', 'TGC':'C', 'TGG':'W', 'TGT':'C', 'TGY':'C', 'TRA':'*', 'TTA':'L', 
	'TTC':'F', 'TTG':'L', 'TTR':'L', 'TTT':'F', 'TTY':'F', 'YTA':'L', 'YTG':'L', 'YTR':'L',
	'---': '-', '...': '-', '~~~': '-'
}


import string
dir(string)

def read_fas(fas_file):
	with open(fas_file) as f:
		seq1 = f.read().split('>')[1:]
		seq2 = [seq.partition('\n') for seq in seq1]
		seq3 = [{seq[0]:seq[2].replace('\n','')} for seq in seq2]
		return seq3
		
#fas=input('your fasta file name: ')
#print(read_fas(fas))

def cds2aa(fas_file,codon = standardWithAmbig):
	all_pro = []
	for seq in read_fas(fas_file):
		k = ''.join(seq.keys())
		v = ''.join(seq.values())
		pro ={}
		if 'U' or 'u' in seq.values():
			v1 = v.translate(v.maketrans('uU','tT'))
			pro[k]=''.join([codon.get(v1.upper()[n:n+3],'X') for n in range(0,len(v1),3)])
		else:
			pro[k]=''.join([codon.get(v.upper()[n:n+3],'X') for n in range(0,len(v),3)])
		all_pro.append(pro)
	return all_pro
	
#print(cds2aa(fas))

def rm_stop_codon(fas_file, codon = standardWithAmbig):
	all_rm = []
	for seq in read_fas(fas_file):
		rm ={}
		k = ''.join(seq.keys())
		v = ''.join(seq.values())
		v1 = v.translate(v.maketrans('uU','tT'))
		for n in range(0,len(v1),3):
			if codon.get(v1.upper()[n:n+3],'X') == '*':
#				print(n)
				rm[k] = v1[0:n]+'---'+v1[n+3:]
				v1 = rm[k]
			else:
				rm[k] = v1
		all_rm.append(rm)
	return all_rm
				
			
#print(rm_stop_codon(fas))	
import glob
for f in glob.glob('*.fasta'):
	op_name = f.replace('fasta','rm_stop_codon.fas')
	with open(op_name,'w') as op:
		dict_list = rm_stop_codon(f)
		for seq in dict_list:
			k = ''.join(seq.keys())
			v = ''.join(seq.values())
			op.write('>'+k+'\n'+v+'\n')
	op_aa_name = f.replace('fasta','aa.fas')
	with open(op_aa_name,'w') as op_aa:
		aa_dict_list = cds2aa(f)
		for aa in aa_dict_list:
			ka = ''.join(aa.keys())
			va = ''.join(aa.values())
			op_aa.write('>'+ka+'\n'+va+'\n') 
	