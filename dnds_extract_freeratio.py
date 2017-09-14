import re
import glob
import pandas as pd
from pandas import DataFrame
"""
n = input("please input the number of species: ")
#dnds_df = open('dnds_free_ratio_df.txt','w')
#dnds_lists_all_terminals = []
dnds_lists_all_terminals = {}
dn_lists_all_terminals = {}
ds_lists_all_terminals = {}
for i in range(n):
    species = 'ctx0000'[:7-len(str(i))]+str(i+1)
#    dnds_lists_all_terminals.append([])
    dnds_lists_all_terminals[species]=[]
    dn_lists_all_terminals[species]=[]
    ds_lists_all_terminals[species]=[]
for op in glob.glob('*codeml'):
#    gene = op.replace(".fa.merge.ctx_codonalign.branch.codeml","")
    gene = re.match(r'.*(ENS\w{12}).*',op).group(1)
    with open(op) as f:
        txt = f.read()
        for i in range(n):
            species = 'ctx0000'[:7-len(str(i))]+str(i+1)
            species_dnds = re.search(r"^w\sratios\sas\slabels\sfor\sTreeView:\n{1}.+"+species+"\s#(\d+\.\d+)",txt,flags = re.M).group(1)
#            species_dnds = re.search(r"w\sratios\sas\slabels\sfor\sTreeView:.+"+species+"\s#(\d+\.\d+)",txt,flags = re.DOTALL|re.S).group(1)
#            dnds_lists_all_terminals[i].append({'id':gene,species:species_dnds})
            dnds_lists_all_terminals[species].append({'id':gene,species:species_dnds})

            species_dn = re.search(r'^dN\stree:\n{1}.+'+species+':\s+(\d+\.\d+)',txt,flags=re.M).group(1)
            dn_lists_all_terminals[species].append({'id':gene,species:species_dn})

            species_ds = re.search(r'^dS\stree:\n{1}.+'+species+':\s+(\d+\.\d+)',txt,flags=re.M).group(1)
            ds_lists_all_terminals[species].append({'id':gene,species:species_ds})

list_contain_sp_dnds_dfs = [DataFrame(sp_dict_lst).set_index('id') for sp_dict_lst in dnds_lists_all_terminals.values()]
list_contain_sp_dn_dfs = [DataFrame(sp_dict_lst).set_index('id') for sp_dict_lst in dn_lists_all_terminals.values()]
list_contain_sp_ds_dfs = [DataFrame(sp_dict_lst).set_index('id') for sp_dict_lst in ds_lists_all_terminals.values()]

dnds_df = pd.concat(list_contain_sp_dnds_dfs,axis=1,ignore_index=False)
dn_df = pd.concat(list_contain_sp_dn_dfs,axis=1,ignore_index=False)
ds_df = pd.concat(list_contain_sp_ds_dfs,axis=1,ignore_index=False)

dnds_df.to_csv('dnds_free_ratio_df_terminal.txt',header=True,index=True)
dn_df.to_csv('dn_free_ratio_df_terminal.txt',header=True,index=True)
ds_df.to_csv('ds_free_ratio_df_terminal.txt',header=True,index=True)
"""

"""
dnds_lists_all_terminals
Out[497]: 
{'ctx0001': [{'ctx0001': '0.3529', 'id': 'ENST00000000412'},
  {'ctx0001': '0.5529', 'id': 'ENST00000000612'}],
 'ctx0002': [{'ctx0002': '0.7764', 'id': 'ENST00000000412'},
  {'ctx0002': '0.8764', 'id': 'ENST00000000612'}],
 'ctx0003': [{'ctx0003': '0.1907', 'id': 'ENST00000000412'},
  {'ctx0003': '0.1907', 'id': 'ENST00000000612'}],
 'ctx0004': [{'ctx0004': '0.0001', 'id': 'ENST00000000412'},
  {'ctx0004': '0.0001', 'id': 'ENST00000000612'}],
 'ctx0005': [{'ctx0005': '0.4190', 'id': 'ENST00000000412'},
  {'ctx0005': '0.4190', 'id': 'ENST00000000612'}],
 'ctx0006': [{'ctx0006': '1.2421', 'id': 'ENST00000000412'},
  {'ctx0006': '1.2421', 'id': 'ENST00000000612'}]}

dnds_lists_all_terminals.values()
Out[498]: 
[[{'ctx0001': '0.3529', 'id': 'ENST00000000412'},
  {'ctx0001': '0.5529', 'id': 'ENST00000000612'}],
 [{'ctx0002': '0.7764', 'id': 'ENST00000000412'},
  {'ctx0002': '0.8764', 'id': 'ENST00000000612'}],
 [{'ctx0003': '0.1907', 'id': 'ENST00000000412'},
  {'ctx0003': '0.1907', 'id': 'ENST00000000612'}],
 [{'ctx0004': '0.0001', 'id': 'ENST00000000412'},
  {'ctx0004': '0.0001', 'id': 'ENST00000000612'}],
 [{'ctx0005': '0.4190', 'id': 'ENST00000000412'},
  {'ctx0005': '0.4190', 'id': 'ENST00000000612'}],
 [{'ctx0006': '1.2421', 'id': 'ENST00000000412'},
  {'ctx0006': '1.2421', 'id': 'ENST00000000612'}]]


list_contain_sp_dfs
Out[500]: 
[                ctx0001
 id                     
 ENST00000000412  0.3529
 ENST00000000612  0.5529,                 ctx0002
 id                     
 ENST00000000412  0.7764
 ENST00000000612  0.8764,                 ctx0003
 id                     
 ENST00000000412  0.1907
 ENST00000000612  0.1907,                 ctx0004
 id                     
 ENST00000000412  0.0001
 ENST00000000612  0.0001,                 ctx0005
 id                     
 ENST00000000412  0.4190
 ENST00000000612  0.4190,                 ctx0006
 id                     
 ENST00000000412  1.2421
 ENST00000000612  1.2421]

dnds_df
Out[501]: 
                ctx0001 ctx0002 ctx0003 ctx0004 ctx0005 ctx0006
id                                                             
ENST00000000412  0.3529  0.7764  0.1907  0.0001  0.4190  1.2421
ENST00000000612  0.5529  0.8764  0.1907  0.0001  0.4190  1.2421
"""


"""
for internal node (ancestor)
"""
"""
dnds_lists_all_internal = {}
dnds_lists_all_internal['ctx0004_ctx0006_ancestor'] = []
#dnds_lists_all_internal['ctx0005_ctx0006_ancestor'] = []

for op in glob.glob('*codeml'):
    geneid = re.match(r'.*(ENS\w{12}).*',op).group(1)
    content = open(op).read()
    dnds1 = re.search(r'^w\sratios\sas\slabels\sfor\sTreeView:\n{1}.+ctx0006\s#\d+\.\d+\s\)\s#(\d+\.\d+)',content,re.M).group(1)
    dnds_lists_all_internal['ctx0004_ctx0006_ancestor'].append({'id':geneid,'ctx0004_ctx0006_ancestor':dnds1})
list_contain_ancestor_dnds_dfs = [DataFrame(ancestor_dict_lst).set_index('id') for ancestor_dict_lst in dnds_lists_all_internal.values()]

dnds_df_internal = pd.concat(list_contain_ancestor_dnds_dfs,axis=1,ignore_index=False)
dnds_df_terminal_plus_internal = pd.concat([dnds_df,dnds_df_internal],axis=1,ignore_index=False)
dnds_df_terminal_plus_internal.to_csv("dnds_df_terminal_plus_internal.csv",header=True,index=True)
"""
"""
dnds_lists_all_internal
Out[557]: 
{'ctx0004_ctx0006_ancestor': [{'ctx0004_ctx0006_ancestor': '3.5393',
   'id': 'ENST00000000412'},
  {'ctx0004_ctx0006_ancestor': '1.5393', 'id': 'ENST00000000612'}]}

list_contain_ancestor_dnds_dfs
Out[558]: 
[                ctx0004_ctx0006_ancestor
 id                                      
 ENST00000000412                   3.5393
 ENST00000000612                   1.5393]
"""


n = input("please input the number of species: ")
dnds_dicts_all_terminals = {}
dn_dicts_all_terminals = {}
ds_dicts_all_terminals = {}
for i in range(n):
    species = 'ctx0000'[:7-len(str(i))]+str(i+1)
    dnds_dicts_all_terminals[species]={}
    dn_dicts_all_terminals[species]={}
    ds_dicts_all_terminals[species]={}
for op in glob.glob('*codeml'):
#    gene = op.replace(".fa.merge.ctx_codonalign.branch.codeml","")
    gene = re.match(r'.*(ENS\w{12}).*',op).group(1)
    with open(op) as f:
        txt = f.read()
        for i in range(n):
            species = 'ctx0000'[:7-len(str(i))]+str(i+1)
            dnds = re.search(r"^w\sratios\sas\slabels\sfor\sTreeView:\n{1}.+"+species+"\s#(\d+\.\d+)",txt,flags = re.M).group(1)
#            species_dnds = re.search(r"w\sratios\sas\slabels\sfor\sTreeView:.+"+species+"\s#(\d+\.\d+)",txt,flags = re.DOTALL|re.S).group(1)
            dnds_dicts_all_terminals[species][gene]=dnds

            dn = re.search(r'^dN\stree:\n{1}.+'+species+':\s+(\d+\.\d+)',txt,flags=re.M).group(1)
            dn_dicts_all_terminals[species][gene]=dn

            ds = re.search(r'^dS\stree:\n{1}.+'+species+':\s+(\d+\.\d+)',txt,flags=re.M).group(1)
            ds_dicts_all_terminals[species][gene]=ds

dnds_df = DataFrame(dnds_dicts_all_terminals)
dn_df = DataFrame(dn_dicts_all_terminals)
ds_df = DataFrame(ds_dicts_all_terminals)

dnds_df.index.name = 'gene_id'
dn_df.index.name = 'gene_id'
ds_df.index.name = 'gene_id'


dnds_df.to_csv('dnds_free_ratio_df_terminal.txt',header=True,index=True)
dn_df.to_csv('dn_free_ratio_df_terminal.txt',header=True,index=True)
ds_df.to_csv('ds_free_ratio_df_terminal.txt',header=True,index=True)



###for internal###
dnds_dicts_all_internal = {}
dnds_dicts_all_internal['ctx0004_ctx0006_ancestor'] = {}
#dnds_dicts_all_internal['ctx0005_ctx0006_ancestor'] = {}

for op in glob.glob('*codeml'):
    geneid = re.match(r'.*(ENS\w{12}).*',op).group(1)
    content = open(op).read()
    dnds1 = re.search(r'^w\sratios\sas\slabels\sfor\sTreeView:\n{1}.+ctx0006\s#\d+\.\d+\s\)\s#(\d+\.\d+)',content,re.M).group(1)
    dnds_dicts_all_internal['ctx0004_ctx0006_ancestor'][geneid] = dnds1

dnds_df_internal = DataFrame(dnds_dicts_all_internal)
dnds_df_internal.index.name = 'gene_id'
dnds_df_terminal_plus_internal = pd.concat([dnds_df,dnds_df_internal],axis=1,ignore_index=False)
dnds_df_terminal_plus_internal.to_csv("dnds_df_terminal_plus_internal.csv",header=True,index=True)

