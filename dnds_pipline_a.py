import pandas as pd
import numpy as np
from pandas import DataFrame, Series

dnds = pd.read_csv("dnds_island_snake.csv",header = 0)
print(dnds.head())
#dnds = dnds.rename(columns = {'Unnamed: 0':'new_name'})
#dnds.columns ['new_names'.....]
#print(dnds.head())
dnds['Unnamed: 0'] = dnds['Unnamed: 0'].replace(['mENST','\.codeml'],['ENST',''],regex=True)
dnds = dnds.rename(columns = {'Unnamed: 0':'Ensembl Transcript ID'})
#print(dnds.head())
mart = pd.read_table('human_mart_export_simp_sort.txt',header = 'infer')
#print(mart.head())
rate_with_go = pd.merge(mart,dnds,how='inner',sort=True)
#print(rate_with_go.head())
print(len(rate_with_go.index))
#102472
rate_with_go_rmdup = rate_with_go.drop_duplicates(keep='first')
print(len(rate_with_go_rmdup.index))
#90842
#df[df.duplicated(['ID'], keep=False)]
dup_rows = rate_with_go[rate_with_go.duplicated(keep='first')]
#print(dup_rows.head())
#print(rate_with_go.duplicated(keep='first').head())
print(len(dup_rows.index))
#11630
GO_uniq = rate_with_go.drop_duplicates(['GO Term Accession'],keep='first')
print(len(GO_uniq.index))
#8500
rate_with_go_rmdup.to_csv('ratebyGO_py.csv',header=True, index=False,na_rep='NA')

GO_mean_rate = rate_with_go_rmdup.groupby('GO Term Accession').mean()
#GO_mean_rate1 = rate_with_go_rmdup.groupby('GO Term Accession')[['Snake_island_viper_6','shedao_heimei_ancestor_11','King_cobra_5']].mean()
#print(GO_mean_rate1.head())
#print(rate_with_go_rmdup.ix[[1,3],2:8].head())
#print(GO_mean_rate.head())
print(len(GO_mean_rate))
#8499
GO_mean_rate.to_csv('GOmeanRate_py.txt',header=True,index=True,sep='\t',na_rep='NA')


import scipy.stats
#go_u_pvalue_output = open('GO_fast_wilcox.txt','w')
#go_u_pvalue_output_sign = open('GO_fast_wilcox_sign.txt','w')
#go_u_pvalue_output_sign_less = open('GO_fast_wilcox_sign_less.txt','w')
#go_u_pvalue_output_less = open('GO_fast_wilcox_less.txt','w')
species_focus = 'Snake_island_viper_6'
species_backgroud = 'shedao_heimei_ancestor_11'
go_u_pvalue = []
go_u_pvalue_less = []
go_u_pvalue_sign = []
go_u_pvalue_sign_less = []
rate_with_go_rmdup_rmna = rate_with_go_rmdup.dropna(subset = [species_focus,species_backgroud])
#len(rate_with_go_rmdup_rmna.index)
#10314
for go_term in rate_with_go_rmdup_rmna['GO Term Accession'].drop_duplicates(keep='first'):
    list1 = rate_with_go_rmdup_rmna[rate_with_go_rmdup_rmna['GO Term Accession'] == go_term][species_focus]
    list2 = rate_with_go_rmdup_rmna[rate_with_go_rmdup_rmna['GO Term Accession'] == go_term][species_backgroud]
#paired    u, pvalue = scipy.stats.wilcoxon(list1,list2)
    u, pvalue = scipy.stats.mannwhitneyu(list1,list2,alternative='greater')
    go_u_pvalue.append([go_term,u,pvalue])
#2342
    if pvalue < 0.05:
        go_u_pvalue_sign.append([go_term,u,pvalue])
    ul, pvaluel = scipy.stats.mannwhitneyu(list1,list2,alternative='less')
    go_u_pvalue_less.append([go_term,ul,pvaluel])
    if pvaluel < 0.05:
        go_u_pvalue_sign_less.append([go_term,ul,pvaluel])
    
#168
#go_u_pvalue_output.write('\n'.join('\t'.join(map(str,go)) for go in go_u_pvalue))
#go_u_pvalue_output_sign.write('\n'.join('\t'.join(map(str,go)) for go in go_u_pvalue_sign))    
#go_u_pvalue_output_sign_less.write('\n'.join('\t'.join(map(str,go)) for go in go_u_pvalue_sign_less))
#go_u_pvalue_output_less.write('\n'.join('\t'.join(map(str,go)) for go in go_u_pvalue_less))

DataFrame(go_u_pvalue).to_csv("GO_fast_wilcox.csv",header=False,index=False,sep='\t')
DataFrame(go_u_pvalue_less).to_csv("GO_fast_wilcox_less.csv",header=False,index=False,sep='\t')
go_u_pvalue_sign_df = DataFrame(go_u_pvalue_sign)
go_u_pvalue_sign_less_df = DataFrame(go_u_pvalue_sign_less)
go_u_pvalue_sign_df.to_csv("GO_fast_wilcox_sign.csv",header=False,index=False,sep='\t')
go_u_pvalue_sign_less_df.to_csv("GO_fast_wilcox_sign_less.csv",header=False,index=False,sep='\t')

"""
intersection bwtween wilcox of python and R
comm -12 <(cut -f1 GO_fast_wilcox_sign.txt |sort |uniq ) <(cut -f1 wilcox_islandviper_vs_shedaoheimeiancestor_signfastGOlabel_Bidirectional.txt|sort |uniq ) 
#146 intersection
difference and intersection 
comm <(cut -f1 GO_fast_wilcox_sign.txt |sort |uniq ) <(cut -f1 wilcox_islandviper_vs_shedaoheimeiancestor_signfastGOlabel_Bidirectional.txt|sort |uniq ) 
"""

"""
rate_with_go_rmdup_rmna = rate_with_go_rmdup.dropna(subset = ['Snake_island_viper_6','Black_brow_viper_4'])
#len(rate_with_go_rmdup_rmna.index)
#10314
for go_term in rate_with_go_rmdup_rmna['GO Term Accession'].drop_duplicates(keep='first'):
    list1 = rate_with_go_rmdup_rmna[rate_with_go_rmdup_rmna['GO Term Accession'] == go_term]['Snake_island_viper_6']
    list2 = rate_with_go_rmdup_rmna[rate_with_go_rmdup_rmna['GO Term Accession'] == go_term]['Black_brow_viper_4']
#paired    u, pvalue = scipy.stats.wilcoxon(list1,list2)
    u, pvalue = scipy.stats.mannwhitneyu(list1,list2,alternative='greater')
    go_u_pvalue.append([go_term,u,pvalue])
#4620
    if pvalue < 0.05:
        go_u_pvalue_sign.append([go_term,u,pvalue])
#7
"""
#go_u_pvalue_sign_df['lable'] = Series(['greater']*len(go_u_pvalue_sign_df), index=go_u_pvalue_sign_df.index)
go_u_pvalue_sign_df.columns = ['GO Term Accession','Mann-Whitney U statistic','p-value']
go_u_pvalue_sign_df_label = go_u_pvalue_sign_df.assign(label = Series(['greater']*len(go_u_pvalue_sign_df)))
go_u_pvalue_sign_less_df.columns = ['GO Term Accession','Mann-Whitney U statistic','p-value']
go_u_pvalue_sign_less_df_label = go_u_pvalue_sign_less_df.assign(label = Series(['less']*len(go_u_pvalue_sign_less_df)))
go_u_pvalue_sign_gl_df = pd.concat([go_u_pvalue_sign_df_label,go_u_pvalue_sign_less_df_label],axis=0,ignore_index=True)
#go_u_pvalue_sign_gl_df_sort = go_u_pvalue_sign_gl_df.sort_values(by='GO Term Accession',ascending=True)

"""
go_u_pvalue_sign_gl_df = go_u_pvalue_sign_df.append(go_u_pvalue_sign_less_df,ignore_index=True)
go_u_pvalue_sign_gl_df = pd.concat([go_u_pvalue_sign_df,go_u_pvalue_sign_less_df],axis=0,ignore_index=True)
go_u_pvalue_sign_gl_df.columns = ['GO Term Accession','Mann-Whitney U statistic','p-value']
go_u_pvalue_sign_gl_df.assign(lable = Series(['greater']*len(go_u_pvalue_sign_df)+['less']*len(go_u_pvalue_sign_less_df)))
"""

go_rate_with_sign_p = pd.merge(GO_mean_rate.reset_index().ix[:,["GO Term Accession",species_focus,species_backgroud]],go_u_pvalue_sign_gl_df,on="GO Term Accession",how='inner',sort=True)
mart_go_uniq_ann = pd.read_table('human_mart_export_GOuniq_sort.txt',header = 'infer',skiprows=[0],skip_blank_lines=True)
go_rate_with_sign_p_ann = pd.merge(mart_go_uniq_ann,go_rate_with_sign_p,on="GO Term Accession",how='inner')

go_rate_with_sign_p_ann_setindex = go_rate_with_sign_p_ann.set_index("GO Term Accession")
active_df_g = go_rate_with_sign_p_ann_setindex[go_rate_with_sign_p_ann_setindex['label']=='greater']
active_df_l = go_rate_with_sign_p_ann_setindex[go_rate_with_sign_p_ann_setindex['label']=='less']
go_rate_with_sign_p_ann_setindex_g = active_df_g[active_df_g[species_focus] > active_df_g[species_backgroud]]
go_rate_with_sign_p_ann_setindex_l = active_df_l[active_df_l[species_focus] < active_df_l[species_backgroud]]
go_rate_with_sign_p_ann_setindex_gl = pd.concat([go_rate_with_sign_p_ann_setindex_g,go_rate_with_sign_p_ann_setindex_l],axis=0,ignore_index=False)
go_rate_with_sign_p_ann_setindex_gl.to_csv("go_rate_with_sign_p_ann_setindex_gl.csv",header=True,index=True)

active_df_true_g_index = active_df_g.index[active_df_g[species_focus] > active_df_g[species_backgroud]]
active_df_true_l_index = active_df_l.index[active_df_l[species_focus] < active_df_l[species_backgroud]]
go_rate_with_sign_p_ann_setindex_gl1 = go_rate_with_sign_p_ann_setindex.ix[active_df_true_g_index|active_df_true_l_index]
#go_rate_with_sign_p_ann_setindex_gl1.to_csv("go_rate_with_sign_p_ann_setindex_gl1.csv",header=True,index=True)

(go_rate_with_sign_p_ann_setindex_gl.sort_index() != go_rate_with_sign_p_ann_setindex_gl1).any(1)
np.where(go_rate_with_sign_p_ann_setindex_gl.sort_index() != go_rate_with_sign_p_ann_setindex_gl1)

"""
pd.DataFrame({'from':df1.values[difference_locations],'to':df2.values[difference_locations]},index = changed.index)
https://stackoverflow.com/questions/17095101/outputting-difference-in-two-pandas-dataframes-side-by-side-highlighting-the-d
In [21]: ne = (df1 != df2).any(1)

In [22]: ne
Out[22]:
0    False
1     True
2     True
dtype: bool
Then we can see which entries have changed:

In [23]: ne_stacked = (df1 != df2).stack()

In [24]: changed = ne_stacked[ne_stacked]

In [25]: changed.index.names = ['id', 'col']

In [26]: changed
Out[26]:
id  col
1   score         True
2   isEnrolled    True
    Comment       True
dtype: bool
Here the first entry is the index and the second the columns which has been changed.

In [27]: difference_locations = np.where(df1 != df2)

In [28]: changed_from = df1.values[difference_locations]

In [29]: changed_to = df2.values[difference_locations]

In [30]: pd.DataFrame({'from': changed_from, 'to': changed_to}, index=changed.index)
Out[30]:
               from           to
id col
1  score       1.11         1.21
2  isEnrolled  True        False
   Comment     None  On vacation
* Note: it's important that df1 and df2 share the same index here. To overcome this ambiguity, you can ensure you only look at the shared labels using df1.index & df2.index, but I think I'll leave that as an exercise.
"""