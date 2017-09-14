import pandas as pd
import glob
import sys
import re
from pandas import DataFrame



if sys.version_info[0] > 2:
    FG_id = input("Given the interesting taxon IDs you want to compare, should be separated by space >")
    BG_id = input("Given the background taxon IDs or ignore it by press enter to use all rests as background >")
    ext = input("Given the input file extension >")
else:
    FG_id = raw_input("Given the interesting taxon IDs you want to compare, should be separated by space >")
    BG_id = raw_input("Given the background taxon IDs or ignore it by press enter to use all rests as background >")
    ext = raw_input("Given the input file extension >")

FG_list = FG_id.split()
BG_list_user = BG_id.split()
"""
FG_id = 'ctx0001 ctx0002 ctx0009'
FG_list = ['ctx0001', 'ctx0002', 'ctx0009']
ext = 'fasta'
"""
ignore_character = ['-','x','X','U','u','O','o','B','b','Z','z','J','j']
AA_list = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V','a','r','n','d','c','e','q','g','h','i','l','k','m','f','p','s','t','w','y','v']
identity_list = []
for file in glob.glob('./*'+ext):
    stander_ID = re.match(r'\S*(ENS\w{12})\S*', file)
    if stander_ID is not None:
        geneid = stander_ID.group(1)
    else:
        geneid = file
    with open(file) as f:
        records = f.read()
        seq1 = records.split('>')[1:]
        seq2 = [seq.partition('\n') for seq in seq1]
        seq3 = {seq[0].strip():seq[2].strip() for seq in seq2}
        all_taxon_list = list(seq3.keys())
        if len(BG_list_user) != 0:
            BG_list = BG_list_user
        else:
            BG_list = list(set(all_taxon_list).difference(set(FG_list)))
        seq_length = len(list(seq3.values())[0])
#        identity_sites = []
        for i in range(0,seq_length):
            FG_sites = []
            BG_sites = []
            for spf in FG_list:
#                if seq3[spf][i] != '-':
                FG_sites.append(seq3[spf][i])
            for spb in BG_list:
#                if seq3[spb][i] != '-':
                BG_sites.append(seq3[spb][i])
            if set(FG_sites).issubset(AA_list) and set(BG_sites) < set(AA_list) and len(set(FG_sites)) == 1 and len(set(FG_sites) & set(BG_sites)) == 0:
                identity_dict = {}
                identity_dict['Gene_ID'] = geneid
                identity_dict['FG_Taxons'] = FG_list
                identity_dict['BG_Taxons'] = BG_list
                identity_dict['FG_AA'] = FG_sites
                identity_dict['Position'] = i+1
                identity_dict['BG_AA'] = BG_sites
                if len(identity_dict) != 0:
                    identity_list.append(identity_dict)

if len(identity_list) != 0:
    identity_df = DataFrame(identity_list)
#    if not identity_df.empty:
    identity_df.set_index('Gene_ID').sort_index().to_csv('_'.join(FG_list)+"identifity_info.csv",header=True,index=True)
else:
    print("No identify sites were found at interesting taxons, you may want to try other taxons")


print(identity_df)





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







