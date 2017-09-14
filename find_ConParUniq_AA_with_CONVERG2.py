import pandas as pd
import glob
import sys
import re
import os
import copy
import subprocess
from pandas import DataFrame



if sys.version_info[0] > 2:
    FG_id = input("Given the two interesting extant taxon IDs you want to compare, should be separated by space >")
    BG_id = input("Given the background taxon IDs or ignore it by press enter to use all rests terminals as background >")
    ANC_id = input("Given the two ancestor IDs corresponding to focused two extant taxons >")
    ext = input("Given the input file extension >")
    prob_test = input("Do you want to compute the probabilities that the observed convergent and parallel substitutions are attributable to random chance? >")
    if prob_test in ['YES', 'yes', 'Y', 'y']:
        model = input("choose the model: 1 for Poisson model and 2 for JTT model >")
else:
    FG_id = raw_input("Given the two interesting extant taxon IDs you want to compare, should be separated by space >")
    BG_id = raw_input("Given the background taxon IDs or ignore it by press enter to use all rests terminals as background >")
    ANC_id = raw_input("Given the two ancestor IDs corresponding to focused two extant taxons >")
    ext = raw_input("Given the input file extension >")
    prob_test = raw_input("Do you want to compute the probabilities that the observed convergent and parallel substitutions are attributable to random chance? >")
    if prob_test in ['YES', 'yes', 'Y', 'y']:
        model = raw_input("choose the model: 1 for Poisson model and 2 for JTT model >")

FG_list = FG_id.split()
BG_list_user = BG_id.split()
ANC_list = ANC_id.split()
ANC_NO = [x.replace("node","") for x in ANC_list]

"""
FG_id = 'ctx0001 ctx0002 ctx0009'
FG_list = ['ctx0001', 'ctx0002', 'ctx0009']
ext = 'fasta'
"""
ignore_character = ['-','x','X','U','u','O','o','B','b','Z','z','J','j']
AA_list = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V','a','r','n','d','c','e','q','g','h','i','l','k','m','f','p','s','t','w','y','v']
unique_share_all = []
convergent_all = []
parallel_all = []

for file in glob.glob('./*'+ext):
    stander_ID = re.match(r'\S*(ENS\w{12})\S*', file)
    if stander_ID is not None:
        geneid = stander_ID.group(1)
    else:
        geneid = os.path.basename(file)
    with open(file) as f:
        records = f.read()
        if 'f' in ext:
            seq1 = records.split('>')[1:]
            seq2 = [seq.partition('\n') for seq in seq1]
            seq3 = {seq[0].strip():seq[2].strip() for seq in seq2}
        elif 'r' in ext:
            all_seqs = re.search(r"List\sof\sextant\sand\sreconstructed\ssequences\n+\s*\d+\s*\d+\n+(.+)\n+Overall\saccuracy\sof",records, re.S)
            FG_NO = [re.search(r"(\d+)_" + x, records).group(1) for x in FG_list]
            if all_seqs is not None:
                seq1 = all_seqs.group(1).strip().split('\n')
                seq2 = [seq.replace(' #', '').partition(" ") for seq in seq1]
                seq3 = {seq[0]: seq[2].strip().replace(' ', '') for seq in seq2}
                if prob_test in ['YES', 'yes', 'Y', 'y']:
                    CONVERG2_fmt = file.replace(ext, 'AA')
                    match_tree = re.search(r"Ancestral\sreconstruction\sby\sAAML\.\n+.+;\n+(.+);\n.+\n+tree\swith\snode\slabels\sfor\sRod\sPage",records, re.S)
                    if match_tree is not None:
                        with open(CONVERG2_fmt, 'w') as CONVERG2_in:
                            seq4 = copy.deepcopy(seq3)
                            for k in list(seq4.keys()):
                                if 'node' in k:
                                    del seq4[k]
                            CONVERG2_in.write("%s  %s\n" % (len(list(seq4.keys())), len(list(seq4.values())[0])))
                            for d in sorted(seq4.keys()):
                                CONVERG2_in.write("%s\n%s\n" % (d,seq4[d]))
                            tree = match_tree.group(1).strip()
                            CONVERG2_in.write("%s\n" % (tree))
                    else:
                        continue
            else:
                continue
        all_taxon_list = list(seq3.keys())
        if len(BG_list_user) != 0:
            BG_list = BG_list_user
        else:
            BG_list = list(x for x in set(all_taxon_list).difference(set(FG_list)) if 'node' not in x)
        seq_length = len(list(seq3.values())[0])
#        identity_sites = []
        count_CONV = 0
        count_PARA = 0
        unique_share_list = []
        convergent_list = []
        parallel_list = []
        for i in range(0,seq_length):
            FG_sites = []
            BG_sites = []
            ANC_sites = []
            for spf in FG_list:
#                if seq3[spf][i] != '-':
                FG_sites.append(seq3[spf][i])
            for spb in BG_list:
#                if seq3[spb][i] != '-':
                BG_sites.append(seq3[spb][i])
            for spa in ANC_list:
                ANC_sites.append(seq3[spa][i])
            if set(FG_sites).issubset(AA_list) and set(BG_sites) < set(AA_list) and len(set(FG_sites)) == 1 and len(set(FG_sites) & set(BG_sites)) == 0:
                unique_share_dict = {}
                unique_share_dict['Gene_ID'] = geneid
                unique_share_dict['FG_Taxons'] = FG_list
                unique_share_dict['BG_Taxons'] = BG_list
                unique_share_dict['ANC_Taxons'] = ANC_list
                unique_share_dict['FG_AA'] = FG_sites
                unique_share_dict['Position'] = i+1
                unique_share_dict['BG_AA'] = BG_sites
                unique_share_dict['ANC_AA'] = ANC_sites
                if len(unique_share_dict) != 0:
                    unique_share_list.append(unique_share_dict)
            if set(FG_sites).issubset(AA_list) and set(BG_sites) < set(AA_list) and len(set(FG_sites)) == 1 and len(set(FG_sites) & set(BG_sites)) == 0 and len(set(ANC_sites)) == len(set(ANC_list)):
#            if set(FG_sites).issubset(AA_list) and set(BG_sites) < set(AA_list) and len(set(FG_sites)) == 1 and len(set(FG_sites) & set(BG_sites)) == 0 and len(set(ANC_sites)) == 2 and len(set(FG_sites) & set(ANC_sites)) == 0:
                count_CONV += 1
                convergent_dict = {}
                convergent_dict['Gene_ID'] = geneid
                convergent_dict['FG_Taxons'] = FG_list
                convergent_dict['BG_Taxons'] = BG_list
                convergent_dict['ANC_Taxons'] = ANC_list
                convergent_dict['FG_AA'] = FG_sites
                convergent_dict['Position'] = i + 1
                convergent_dict['BG_AA'] = BG_sites
                convergent_dict['ANC_AA'] = ANC_sites
                if len(convergent_dict) != 0:
                    convergent_list.append(convergent_dict)
            if set(FG_sites).issubset(AA_list) and set(BG_sites) < set(AA_list) and len(set(FG_sites)) == 1 and len(set(FG_sites) & set(BG_sites)) == 0 and len(set(ANC_sites)) == 1 and len(set(FG_sites) & set(ANC_sites)) == 0:
                count_PARA += 1
                parallel_dict = {}
                parallel_dict['Gene_ID'] = geneid
                parallel_dict['FG_Taxons'] = FG_list
                parallel_dict['BG_Taxons'] = BG_list
                parallel_dict['ANC_Taxons'] = ANC_list
                parallel_dict['FG_AA'] = FG_sites
                parallel_dict['Position'] = i + 1
                parallel_dict['BG_AA'] = BG_sites
                parallel_dict['ANC_AA'] = ANC_sites
                if len(parallel_dict) != 0:
                    parallel_list.append(parallel_dict)

        if prob_test in ['YES', 'yes', 'Y', 'y']:
            if count_PARA != 0 or count_CONV != 0:
                p = subprocess.Popen(["converg2", CONVERG2_fmt], shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
                with p.stdin:
                    for v in [str(model) + "\n", ANC_NO[0] + "\n", FG_NO[0] + "\n", ANC_NO[1] + "\n", FG_NO[1] + "\n", str(count_PARA) + "\n", str(count_CONV) + "\n"]:
                        p.stdin.write(v)
                p.stdin.close()
                CONVERG2_out = p.stdout.read()
                exp_para = re.search(r"Expected\snumber\sof\sparallel-change\ssites\s*=\s*(\d*\.?\d*)", CONVERG2_out,re.S).group(1)
                obs_para = re.search(r"Observed\snumber\sof\sparallel-change\ssites\s*=\s*(\d*\.?\d*)", CONVERG2_out,re.S).group(1)
                prob_para = re.findall(r"probability\s*=\s*(\d*\.?\d*)", CONVERG2_out, re.I)[0]
                exp_conv = re.search(r"Expected\snumber\sof\sconvergent-change\ssites\s*=\s*(\d*\.?\d*)", CONVERG2_out,re.S).group(1)
                obs_conv = re.search(r"Observed\snumber\sof\sconvergent-change\ssites\s*=\s*(\d*\.?\d*)", CONVERG2_out,re.S).group(1)
                prob_conv = re.findall(r"probability\s*=\s*(\d*\.?\d*)", CONVERG2_out, re.I)[1]
                for dc in convergent_list:
                    dc["Exp_CONV"] = exp_conv
                    dc["Obs_CONV"] = obs_conv
                    dc["Prob_CONV"] = prob_conv
                for dp in parallel_list:
                    dp["Exp_PARA"] = exp_para
                    dp["Obs_PARA"] = obs_para
                    dp["Prob_PARA"] = prob_para
        unique_share_all.extend(unique_share_list)
        convergent_all.extend(convergent_list)
        parallel_all.extend(parallel_list)
if len(unique_share_all) != 0:
    unique_share_df = DataFrame(unique_share_all)
#    if not identity_df.empty:
    unique_share_df.set_index('Gene_ID').sort_index().to_csv('_'.join(FG_list)+"_uniquely_shared_info.csv",header=True,index=True)
else:
    print("No uniquely shared sites were found at interesting taxons, you may want to try other taxons")

if len(convergent_all) != 0:
    convergent_df = DataFrame(convergent_all)
#    if not identity_df.empty:
    convergent_df.set_index('Gene_ID').sort_index().to_csv('_'.join(FG_list)+"_convergent_info_prob.csv",header=True,index=True)
else:
    print("No convergent sites were found at interesting taxons, you may want to try other taxons")

if len(parallel_all) != 0:
    parallel_df = DataFrame(parallel_all)
#    if not identity_df.empty:
    parallel_df.set_index('Gene_ID').sort_index().to_csv('_'.join(FG_list)+"_parallel_info_prob.csv",header=True,index=True)
else:
    print("No parallel sites were found at interesting taxons, you may want to try other taxons")



print(BG_list)
print(ANC_list)
print(FG_list)





"""
https://stackoverflow.com/questions/3844801/check-if-all-elements-in-a-list-are-identical
def checkEqual1(iterator):
    iterator = iter(iterator)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == rest for rest in iterator)
One-liner:

def checkEqual2(iterator):
   return len(set(iterator)) <= 1
Also one-liner:

def checkEqual3(lst):
   return lst[1:] == lst[:-1]

# http://stackoverflow.com/q/3844948/
def checkEqualIvo(lst):
    return not lst or lst.count(lst[0]) == len(lst)

# http://stackoverflow.com/q/3844931/
def checkEqual6502(lst):
    return not lst or [lst[0]]*len(lst) == lst
"""







