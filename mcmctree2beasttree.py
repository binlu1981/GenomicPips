import re
import glob

for file in glob.glob('mcmctree*nwk'):
	op = file.replace('nwk','nex')
	with open(file) as f:
		trees = f.read()
		beast_tree = re.sub(r'(\d*\.?\d*)\-(\d*\.?\d*)',r'[&95%={\1, \2}]',trees)
	with open(op,'w') as o:
		o.write("#NEXUS\nBEGIN TREES;\n\n\tUTREE 1 = "+beast_tree+"END;\n")
	
