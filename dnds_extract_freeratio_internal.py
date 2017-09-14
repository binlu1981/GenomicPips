import re
import glob
import pandas as pd
from pandas import DataFrame


dnds_df=pd.read_csv('dnds_free_ratio_df_terminal.txt',header=0)
#dn_df.to_csv('dn_free_ratio_df_terminal.txt',header=True,index=True)
#ds_df.to_csv('ds_free_ratio_df_terminal.txt',header=True,index=True)



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

