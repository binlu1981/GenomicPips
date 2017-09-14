import os
import shutil

shareid = open("shareID.txt",'w')

"""make a list contain mutiple dicts (key:foldname,value:filename in folders)"""
cwd = os.getcwd()
# print cwd
# print os.path.split(cwd)
# print os.listdir(cwd)
# for root, dirs, files in os.walk(cwd):
# folderlist = [d for d in os.listdir(cwd) if os.path.isdir(os.path.join(cwd, d))]
folderlist = next(os.walk(cwd))[1]
# print folderlist
filediclist = []
for folder in folderlist:
	dic = {}
	filelist = os.listdir(os.path.join(cwd,folder))
	# print os.path.join(cwd,folder)
	dic[folder] = filelist
	filediclist.append(dic)
# print filediclist

"""find the intersection id name"""
filelistslist = []
for d in filediclist:
	# print d.values()	
	# filelistslist+=d.values()
	filelistslist.extend(d.values())
# print filelistslist
# sharefile = set.intersection(*map(set,filelistslist))
sharefile = reduce(set.intersection, map(set, filelistslist))
sharefilelist = list(sharefile)
sharefileid = ("\n").join(sharefilelist)
print sharefileid
shareid.write("%s" % sharefileid)

"""open the shard files in each folder and combine shared files"""
for id in sharefilelist:
	pathlist = []
	with open(id,'w') as combinefile:
		for folder in folderlist:
			pathlist.append(os.path.join(cwd,folder,id))
		for sharefile in pathlist:
			# print file
			tempfile = open(sharefile).read()
			if tempfile[-1] != "\n":
				tempfile+="\n"
			combinefile.write("%s" % tempfile)
		combinefile.close()



