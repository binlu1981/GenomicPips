import re
import glob

for pb in glob.glob("*.prob.txt"):
#	print pb
	focustaxon=pb.replace("prob.txt","")
#	parallelsingnificant=open(focustaxon+"parallelsingnificant.txt","w")
	convergesingnificant=open(focustaxon+"convergesingnificant.txt","w")
	convergeoutput=open(pb).read()
	#print convergeoutput
	eplst=convergeoutput.split("\n\n\n")
#	print eplst
#	for ep in eplst:
	pattern3 = re.compile("(\w+)\n(\w+)\n.+\n\n.+\n.+\n.+\n(Observed\snumber\sof\sparallel-change)\ssites\s\=\s\d+\nprobability=\s(\d+\.\d+)")
	pattern4 = re.compile("(\w+)\n(\w+)\n.+\n\n.+\n.+\n.+\n.+\n.+\n\n.+\n(Observed\snumber\sof\sconvergent-change)\ssites\s\=\s\d+\nprobability=\s(\d+\.\d+)")
	#pattern3 = re.compile("(\w+)\r\n(\w+)\r\n.+\r\n\r\n.+\r\n.+\r\n.+\r\n(Observed\snumber\sof\sparallel-change)\ssites\s\=\s\d+\r\nProbability=\s(\d+\.\d+)")
	#pattern4 = re.compile("(\w+)\r\n(\w+)\r\n.+\r\n\r\n.+\r\n.+\r\n.+\r\n.+\r\n.+\r\n\r\n.+\r\n(Observed\snumber\sof\sconvergent-change)\ssites\s\=\s\d+\r\nprobability=\s(\d+\.\d+)")
	for ep in eplst:
		parallelprob = re.findall(pattern3,ep)
#		print ep
#		print parallelprob
#		print len(parallelprob)
		for p in range(0,len(parallelprob)):
			pvalue=parallelprob[p][3]
			if float(pvalue) < 0.05:
#				print parallelprob[p]
				convergesingnificant.write("%s:\t%s:\t%s:\t%s\n" % (parallelprob[p]))
		convergeprob = re.findall(pattern4,ep)
#		print convergeprob
#		print len(convergeprob)
		for c in range(0,len(convergeprob)):
			pvaluec=convergeprob[c][3]
			if float(pvaluec) < 0.05:
#				print convergeprob[c]
				convergesingnificant.write("%s:\t%s:\t%s:\t%s\n" % (convergeprob[c]))
