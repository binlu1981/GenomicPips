import pandas as pd
import glob
import sys
import re
import os
import shutil
import copy
from collections import Counter
from pandas import DataFrame

def fas2dict(fas_records):
    seq1 = fas_records.split('>')[1:]
    seq2 = [seq.partition('\n') for seq in seq1]
    seq3 = {seq[0].strip(): seq[2].strip() for seq in seq2}
    return seq3

def rst2dict(rst_records):
    seq1 = rst_records.group(1).strip().split('\n')
    seq2 = [seq.replace(' #', '').partition(" ") for seq in seq1]
    seq3 = {seq[0]: seq[2].strip().replace(' ', '') for seq in seq2}
    return seq3

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

def rate2rat(rat_records):
    match_content = re.search(r"(\s+Site\s+Freq\s+Data\s+Rate.+)\n+lnL\s+=", rat_records, re.S)
    lst = match_content.group(1).strip().replace("\n\n", "\n").split('\n')
    lst1 = [line.strip() for line in lst]
    lst2 = [re.sub(r"\s+", ",", line) for line in lst1]
    mat = "\n".join(lst2).replace("posterior,mean,&,category", "posterior_mean&category")
    df = pd.read_csv(StringIO(mat), sep=",", header=0)
    return df

def rst2siteFreq(seqs_dict):
    siteFreq_dict = {}
    taxon = list(seqs_dict.keys())
    seq_len = len(list(seqs_dict.values())[0])
    seq_num = len(taxon)
    for i in range(0, seq_len):
        each_site_aa_list = []
        for sp in sorted(taxon):
            each_site_aa_list.append(seqs_dict[sp][i])
        each_site_count = Counter(each_site_aa_list)
        aaFreq_each_site = {}
        for aa in AA_list[:20]:
            if aa in each_site_count:
                aaFreq_each_site[aa] = each_site_count[aa]/seq_num
            else:
                aaFreq_each_site[aa] = 0
        siteFreq_dict[i] = aaFreq_each_site
    df = DataFrame(siteFreq_dict).T
    return df



def creat_folder(folders_list):
    for x in folders_list:
        folder_path = os.path.join(os.getcwd(),x)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print("Created: "+x)


if sys.version_info[0] > 2:
    ext = input("Given the input file extension >")
    from io import StringIO
else:
    ext = raw_input("Given the input file extension >")
    from StringIO import StringIO

folders = ["01_NodeSeq","02_Trees","03_rate","04_siteFreq"]
creat_folder(folders)
work_paths = [os.path.join(os.getcwd(),x) for x in folders]


ignore_character = ['-','x','X','U','u','O','o','B','b','Z','z','J','j']
AA_list = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V','a','r','n','d','c','e','q','g','h','i','l','k','m','f','p','s','t','w','y','v']

for file in glob.glob('./*'+ext):
#    site_rates = file.replace(ext,'rat')
    NodeSeq_fmt = file.replace(ext,'nodes')
    Tree_fmt = file.replace(ext,'ntree')
    rate_fmt = file.replace(ext, 'rat')
    siteFreq_fmt = file.replace(ext,"freq")
    stander_ID = re.match(r'\S*(ENS\w{12})\S*', file)
    if stander_ID is not None:
        geneid = stander_ID.group(1)
    else:
        geneid = os.path.basename(file)
    with open(file) as f:
        records = f.read()
        if 'f' in ext:
            seqs_dict = fas2dict(records)
        elif 'r' in ext:
            all_seqs = re.search(r"List\sof\sextant\sand\sreconstructed\ssequences\n+\s*\d+\s*\d+\n+(.+)\n+Overall\saccuracy\sof",records, re.S)
            if all_seqs is not None:
                seqs_dict = rst2dict(all_seqs)
            else:
                continue
            with open(os.path.join(work_paths[0],NodeSeq_fmt),"w") as NodeSeq:
                for k in sorted(seqs_dict.keys()):
                    NodeSeq.write("%s\t%s\n" % (k,seqs_dict[k]))
            with open(os.path.join(work_paths[1],Tree_fmt),"w") as Trees:
                Trees.write("%s\n" % (rst2tree(records)))
            siteFreq_df = rst2siteFreq(seqs_dict)
            siteFreq_df.to_csv(os.path.join(work_paths[3],siteFreq_fmt),index=False,header=False,sep="\t",float_format='%.16f')
    with open(rate_fmt) as rt:
        rat_record = rt.read()
        rate_df = rate2rat(rat_record)
        rate_df["Rate"].to_csv(os.path.join(work_paths[2],rate_fmt),index=False,float_format='%.3f')











