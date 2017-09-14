import glob
import re
files = glob.glob('*codeml')
output = open('not_saturated_genes.txt','w')
not_saturated_list = []
for f in files:
	content = open(f).read()
#	for line in content:
	seq_len = float(re.match('\s+\d+\s+(\d+)\n',content).group(1))
	seq_len_P = float(re.search('\s+\d+\s+(\d+)\s+P\n',content).group(1))
	seq_n = float(re.match('\s+(\d)+\s+\d+\n',content).group(1))
	N = float(re.search('S\*dS\n\n\s+\d+\.\.\d+\s+\d+\.\d+\s+(\d+\.\d+)\s+(\d+\.\d+)\s+',content,re.S|re.M).group(1))
#	pattern = re.compile('S\*dS\n\n\s+\d+\.\.\d+\s+\d+\.\d+\s+(\d+\.\d+)\s+(\d+\.\d+)\s+',re.S | re.M)
#	N = pattern.search(content).group)(1)
	S = float(re.search('S\*dS\n\n\s+\d+\.\.\d+\s+\d+\.\d+\s+(\d+\.\d+)\s+(\d+\.\d+)\s+',content,re.S|re.M).group(2))
	dN_avg = float(re.search('tree\s+length\s+for\s+dN\:\s+(\d+\.\d+)',content).group(1))/(2*seq_n-3)
	dS_avg = float(re.search('tree\s+length\s+for\s+dS\:\s+(\d+\.\d+)',content).group(1))/(2*seq_n-3)
	NdN_avg = N*dN_avg
	SdS_avg = S*dS_avg 
#	print(seq_len,seq_len_P,seq_n,N,S,dN_avg,dS_avg,NdN_avg,SdS_avg)
	if dS_avg <= 2.0 and NdN_avg >= 1 and SdS_avg >=1 and seq_len >= 150.0:
		not_saturated_list.append(f.split('.')[0])
output.write('\n'.join(not_saturated_list))
output.close()
		