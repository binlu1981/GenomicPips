import glob
import re
import os
import sys




def rst2tree(rst_records):
    tree_branch = re.search(r"Ancestral\sreconstruction\sby\sAAML\.\n+(.+);\n.+;.+\n+tree\swith\snode\slabels\sfor\sRod\sPage", rst_records,re.S).group(1).strip()
    tree_topology_anc = re.search(r"tree\swith\snode\slabels\sfor\sRod\sPage\'s\sTreeView\n(.+);", rst_records, re.S).group(1).strip()
    anc_nodes = re.findall(r"\)\s(\d+)", tree_topology_anc)
    for anc in anc_nodes:
        temp_tree = re.sub("\):", ")" + anc + ":", tree_branch, 1)
        tree_branch = temp_tree
    tree_branch = tree_branch + anc_nodes[-1] + ":0.00000;"
    tree_branch_anc = re.sub(r"\s+", "", tree_branch)
    return tree_branch_anc




def creat_folder(folders_list):
    for x in folders_list:
        folder_path = os.path.join(os.getcwd(),x)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print("Created: "+x)

if sys.version_info[0] > 2:
    ext = input("Given the input file extension >")
else:
    ext = raw_input("Given the input file extension >")

#folders = ["01_NodeSeq","02_Trees","03_rate","04_siteFreq"]
#creat_folder(folders)
#work_paths = [os.path.join(os.getcwd(),x) for x in folders]

for file in glob.glob('./*'+ext):
    Tree_fmt = file.replace(ext,'ntree')
    stander_ID = re.match(r'\S*(ENS\w{12})\S*', file)
    if stander_ID is not None:
        geneid = stander_ID.group(1)
    else:
        geneid = os.path.basename(file)
    with open(file) as f:
        records = f.read()
        with open(Tree_fmt,"w") as Trees:
            Trees.write("%s\n" % (rst2tree(records)))











