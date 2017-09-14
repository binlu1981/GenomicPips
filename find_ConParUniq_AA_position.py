import pandas as pd
import glob
import sys
import re
import os
from pandas import DataFrame



if sys.version_info[0] > 2:
    FG_id = input("Given the two or more interesting extant taxon IDs you want to compare, should be separated by space >")
    BG_id = input("Given the background taxon IDs or ignore it by press enter to use all rests terminals as background >")
    ANC_id = input("Given the two ancestor IDs corresponding to focused two extant taxons >")
    ext = input("Given the input file extension >")
else:
    FG_id = raw_input("Given the two or more interesting extant taxon IDs you want to compare, should be separated by space >")
    BG_id = raw_input("Given the background taxon IDs or ignore it by press enter to use all rests terminals as background >")
    ANC_id = raw_input("Given the two ancestor IDs corresponding to focused two extant taxons >")
    ext = raw_input("Given the input file extension >")

FG_list = FG_id.split()
BG_list_user = BG_id.split()
ANC_list = ANC_id.split()
"""
FG_id = 'ctx0001 ctx0002 ctx0009'
FG_list = ['ctx0001', 'ctx0002', 'ctx0009']
ext = 'fasta'
"""
ignore_character = ['-','x','X','U','u','O','o','B','b','Z','z','J','j']
AA_list = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V','a','r','n','d','c','e','q','g','h','i','l','k','m','f','p','s','t','w','y','v']
unique_share_list = []
convergent_list = []
parallel_list = []
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
            if all_seqs is not None:
                seq1 = all_seqs.group(1).strip().split('\n')
                seq2 = [seq.replace(' #', '').partition(" ") for seq in seq1]
                seq3 = {seq[0]: seq[2].strip().replace(' ', '') for seq in seq2}
            else:
                continue
        all_taxon_list = list(seq3.keys())
        if len(BG_list_user) != 0:
            BG_list = BG_list_user
        else:
            BG_list = list(x for x in set(all_taxon_list).difference(set(FG_list)) if 'node' not in x)
        seq_length = len(list(seq3.values())[0])
#        identity_sites = []
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

if len(unique_share_list) != 0:
    unique_share_df = DataFrame(unique_share_list)
#    if not identity_df.empty:
    unique_share_df.set_index('Gene_ID').sort_index().to_csv('_'.join(FG_list)+"_uniquely_shared_info.csv",header=True,index=True)
else:
    print("No uniquely shared sites were found at interesting taxons, you may want to try other taxons")

if len(convergent_list) != 0:
    convergent_df = DataFrame(convergent_list)
#    if not identity_df.empty:
    convergent_df.set_index('Gene_ID').sort_index().to_csv('_'.join(FG_list)+"_convergent_info.csv",header=True,index=True)
else:
    print("No convergent sites were found at interesting taxons, you may want to try other taxons")

if len(parallel_list) != 0:
    parallel_df = DataFrame(parallel_list)
#    if not identity_df.empty:
    parallel_df.set_index('Gene_ID').sort_index().to_csv('_'.join(FG_list)+"_parallel_info.csv",header=True,index=True)
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







