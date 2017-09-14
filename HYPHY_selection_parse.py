# -*- coding: utf-8 -*-
import pandas as pd
import re
import glob
from pandas import DataFrame, Series
import numpy as np
import json
import sys
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO


method_dict = {"f":"FUBAR","a":"aBSREL","b":"BUSTED","r":"RELAX","m":"MEME","s":"SLAC","fel":"FEL"}
print(method_dict)
if sys.version_info[0] >= 3:
    method = input("choose your method? ")
    while method not in method_dict:
        method = input("Error! Note that it is case sensitive,retry ")
else:
    method = raw_input("choose your method? ")
#method = 'a'
#json_ext = ["a","b","r"]
    while method not in method_dict:
        method = raw_input("Error! Note that it is case sensitive,retry ")
    

    
def parsing_aBSREL(p_threshold = 0.05,ensembl_id = True,json_and_screen = True):
#    selection_gene_dict = {}
    selection_gene_dict_lst = []    
    for file in glob.glob('*aBSREL*.json'):
        if ensembl_id == True:
            geneid = re.match(r'\S*(ENS\w{12})\S*',file).group(1)
        else:
            geneid = file
        with open(file) as f:
            records = json.loads(f.read())
            result = records["test results"]
            for sp_id in result:
                if (result[sp_id]["p"] < p_threshold):
                    selection_sp_dict = {}
#                    selection_sp_dict[sp_id] = result[sp_id]
#                    selection_gene_dict[geneid] = selection_sp_dict
                    selection_sp_dict["Taxon"] = sp_id
                    selection_sp_dict["LRT-value"] = result[sp_id]["LRT"]
                    selection_sp_dict["p-value"] = result[sp_id]["p"]
                    selection_sp_dict["uncorrected_p-value"] = result[sp_id]["uncorrected p"]
                    selection_sp_dict["Gene_ID"] = geneid
                    selection_gene_dict_lst.append(selection_sp_dict)
    if len(selection_gene_dict_lst):
        selection_df = DataFrame(selection_gene_dict_lst)
    else:
        selection_df = pd.DataFrame(columns=['Gene_ID','values'])    
    if json_and_screen == True:
#        pd.options.display.float_format = '{:,.7f}'.format
#        selection_gene_dict1 = {}
        selection_gene_df_lst = []
        for file1 in glob.glob('*aBSREL*screen*'):
            if ensembl_id == True:
                geneid1 = re.match(r'\S*(ENS\w{12})\S*',file1).group(1)
            else:
                geneid1 = file1
            with open(file1) as f1:
                records1 = f1.read()
#                lst = re.findall(r"(\S*)\sp\s=\s(-?\d+\.?\d*\S*)",records1) #do not use '-' as speices id
#                selection_gene_dict1[geneid1] = [{x[0]:x[1]} for x in lst]
                ctt_str = re.search(r"Summary\sof\sbranches\sunder\sepisodic\sselection.+:\n(.+)\s=\=\=\sCPU\sTIME\sREPORT",records1,flags=re.S).group(1)
#                if (ctt_str != "No branches were found to be under selection at p <= 0.05")&(ctt_str is not None):
                if len(ctt_str):
                    if ctt_str.strip() != "No branches were found to be under selection at p <= 0.05":
                        ctt_str1 = ctt_str.replace("\t","").replace(" p = ",",").strip()
                        df_each_gene = pd.read_csv(StringIO(ctt_str1),sep=',',header=None,names=["Taxon","p-value"])
                        df_each_gene = df_each_gene.assign(Gene_ID = Series(geneid1, index=df_each_gene.index))
                        selection_gene_df_lst.append(df_each_gene)
        if len(selection_gene_df_lst):
            selection_df1 = pd.concat(selection_gene_df_lst,axis=0,ignore_index=True)
        else:
            selection_df1 = pd.DataFrame(columns=['Gene_ID','values'])
#        return {"screen":selection_gene_dict1,"json":selection_gene_dict}
        return {"screen":selection_df1,"json":selection_df}
    else:
        return {"json":selection_df}




def parsing_BUSTED(p_threshold = 0.05,ensembl_id = True,json_and_screen = True):        
    selection_gene_dict_lst_full = []
    selection_gene_dict_lst_foreground = []
    for file in glob.glob('*BUSTED*.json'):
        if ensembl_id == True:
            geneid = re.match(r'\S*(ENS\w{12})\S*',file).group(1)
        else:
            geneid = file
        with open(file) as f:
            records = json.loads(f.read())
            if records["test results"]["p"] < p_threshold:
                if records["fits"]["Unconstrained model"]["rate distributions"].get("BG") is not None:
                    selection_index_dict_foreground = {}
                    selection_index_dict_foreground["Gene_ID"] = geneid
                    selection_index_dict_foreground["FG_omega"] = records["fits"]["Unconstrained model"]["rate distributions"]["FG"][2][0]
                    selection_index_dict_foreground["FG_weight"] = records["fits"]["Unconstrained model"]["rate distributions"]["FG"][2][1]
                    selection_index_dict_foreground["p-value"] = records["test results"]["p"]
                    selection_index_dict_foreground["Unconstrained_log-likelihood"] = records["fits"]["Unconstrained model"]["log-likelihood"]
                    selection_index_dict_foreground["constrained_log-likelihood"] = records["fits"]["Constrained model"]["log-likelihood"]
                    selection_index_dict_foreground["BG_omega"] = records["fits"]["Unconstrained model"]["rate distributions"]["BG"][2][0]
                    selection_index_dict_foreground["BG_weight"] = records["fits"]["Unconstrained model"]["rate distributions"]["BG"][2][1]
                    selection_gene_dict_lst_foreground.append(selection_index_dict_foreground)
                else:
                    selection_index_dict_full = {}
                    selection_index_dict_full["Gene_ID"] = geneid
                    selection_index_dict_full["FG_omega"] = records["fits"]["Unconstrained model"]["rate distributions"]["FG"][2][0]
                    selection_index_dict_full["FG_weight"] = records["fits"]["Unconstrained model"]["rate distributions"]["FG"][2][1]
                    selection_index_dict_full["p-value"] = records["test results"]["p"]
                    selection_index_dict_full["Unconstrained_log-likelihood"] = records["fits"]["Unconstrained model"]["log-likelihood"]
                    selection_index_dict_full["constrained_log-likelihood"] = records["fits"]["Constrained model"]["log-likelihood"]
                    selection_gene_dict_lst_full.append(selection_index_dict_full)
    if len(selection_gene_dict_lst_foreground):
        selection_df_foreground = DataFrame(selection_gene_dict_lst_foreground)
    else:
        selection_df_foreground = pd.DataFrame(columns=['Gene_ID','values'])    
    if len(selection_gene_dict_lst_full):
        selection_df_full = DataFrame(selection_gene_dict_lst_full)
    else:
        selection_df_full = pd.DataFrame(columns=['Gene_ID','values'])        
    if json_and_screen == True:
        selection_gene_dict_lst1 = []
        for file1 in glob.glob('*BUSTED*screen*'):
            if ensembl_id == True:
                geneid1 = re.match(r'\S*(ENS\w{12})\S*',file1).group(1)
            else:
                geneid1 = file1
            with open(file1) as f1:
                records1 = f1.read()
                lst1 = re.findall("Log\(L\)\s\=\s(-?\d+\.?\d*\S*)\.\sUnrestricted\sclass\somega\s\=\s(-?\d+\.?\d*\S*)\s\(weight\s\=\s(-?\d+\.?\d*\S*)\).+Log\(L\)\s\=\s(-?\d+\.?\d*\S*).+p\s\=\s(-?\d+\.?\d*\S*)",records1,flags=re.S)
                if len(lst1) != 0:
                    if float(lst1[0][4]) < p_threshold:
                        selection_index_dict1 ={}
                        selection_index_dict1["Gene_ID"] = geneid1
                        selection_index_dict1["omega"] = lst1[0][1]
                        selection_index_dict1["weight"] = lst1[0][2]
                        selection_index_dict1["p-value"] = lst1[0][4]
                        selection_index_dict1["Unconstrained_log-likelihood"] = lst1[0][0]
                        selection_index_dict1["constrained_log-likelihood"] = lst1[0][3]
                        selection_gene_dict_lst1.append(selection_index_dict1)
        if len(selection_gene_dict_lst1):
            selection_df1 = DataFrame(selection_gene_dict_lst1)
        else:
            selection_df1 = pd.DataFrame(columns=['Gene_ID','values'])    
        return {"screen":selection_df1,"json_foreground":selection_df_foreground,"json_full":selection_df_full}                
    else:
        return {"json_foreground":selection_df_foreground,"json_full":selection_df_full}




    
    
    
    
def parsing_RELAX(p_threshold = 0.05,ensembl_id = True,json_and_screen = True):
    selection_gene_lst = [] 
    for file in glob.glob("*RELAX*.json"):
        if ensembl_id == True:
            geneid = re.match(r'\S*(ENS\w{12})\S*',file).group(1)
        else:
            geneid = file
        with open(file) as f:
            records = json.loads(f.read())
            if records["relaxation-test"]["p"] < p_threshold:
                selection_index_dict = {}
                selection_index_dict["Gene_ID"] = geneid
                selection_index_dict["K_value"] = records["fits"]["Alternative"]["K"]
                selection_index_dict["p-value"] = records["relaxation-test"]["p"]
                selection_index_dict["Null_log-likelihood"] = records["fits"]["Null"]["log-likelihood"]
                selection_index_dict["Alternative_log-likelihood"] = records["fits"]["Alternative"]["log-likelihood"]
                selection_gene_lst.append(selection_index_dict)
    if len(selection_gene_lst):
        selection_df = DataFrame(selection_gene_lst)
    else:
        selection_df = pd.DataFrame(columns=['Gene_ID','values'])    
    if json_and_screen == True:
        selection_gene_lst1 = []
        for file1 in glob.glob("*RELAX*screen*"):            
            if ensembl_id == True:
                geneid1 = re.match(r'\S*(ENS\w{12})\S*',file1).group(1)
            else:
                geneid1 = file1
            with open(file1) as f1:
                records1 = f1.read()
                lst1 = re.findall(r"Fitting\sthe\sRELAX\snull\smodel.+Log\(L\)\s\=\s(-?\d+\.?\d*\S*).+Log\(L\)\s\=\s(-?\d+\.?\d*\S*)\.\sRelaxation\sparameter\sK\s\=\s(-?\d+\.?\d*\S*).+p\s\=\s(-?\d+\.?\d*\S*)",records1,flags=re.S)
                if len(lst1) != 0:
                    if float(lst1[0][3]) < p_threshold:
                        selection_index_dict1 = {}
                        selection_index_dict1["Gene_ID"] = geneid1
                        selection_index_dict1["Null_log-likelihood"] = lst1[0][0]
                        selection_index_dict1["Alternative_log-likelihood"] = lst1[0][1]
                        selection_index_dict1["K_value"] = lst1[0][2]
                        selection_index_dict1["p-value"] = lst1[0][3]                    
                        selection_gene_lst1.append(selection_index_dict1)
        if len(selection_gene_lst1):
            selection_df1 = DataFrame(selection_gene_lst1)
        else:
            selection_df1 = pd.DataFrame(columns=['Gene_ID','values'])    
        return {"screen":selection_df1,"json":selection_df}                
    else:
        return {"json":selection_df}

            
                        
                                    

    
    
    
def parsing_FUBAR(Prob = 0.9,ensembl_id = True,csv_and_screen = True):    
    selection_gene_df_lst = []
    for file in glob.glob("*fubar.csv"):
        if ensembl_id == True:
            geneid = re.match(r'\S*(ENS\w{12})\S*',file).group(1)
        else:
            geneid = file
        with open(file) as f:
            records = pd.read_csv(f,sep=',',header='infer')
#            active_index = records.index[records["Prob[alpha<beta]"]>= Prob]
#            active_data = records.ix[active_index]
            active_data = records[records["Prob[alpha<beta]"]>= Prob]
            active_data = active_data.assign(Gene_ID = Series(geneid, index=active_data.index))
#            active_data["Gene_ID"] = Series(geneid, index=active_data.index)
#            active_data.loc[:,"Gene_ID"] = Series(geneid, index=active_data.index)
            if len(active_data):
                selection_gene_df_lst.append(active_data)
    if len(selection_gene_df_lst):
        selection_dfs = pd.concat(selection_gene_df_lst,axis=0,ignore_index=True)
    else:
        selection_dfs = pd.DataFrame(columns=['Gene_ID','values'])
    if csv_and_screen == True:
        selection_gene_df_lst1 = []
        for file1 in glob.glob("*fubar*screen*"):
#            selection_index_lst={}
            if ensembl_id == True:
                geneid1 = re.match(r'\S*(ENS\w{12})\S*',file1).group(1)
            else:
                geneid1 = file1
            with open(file1) as f1:
                records1 = f1.read()
#                match_content = re.search(r"N_eff\n{1}(\d+)\s(-?\d+\.?\d*\S*)\s(-?\d+\.?\d*\S*)\s(-?\d+\.?\d*\S*)\s(-?\d+\.?\d*\S*)",records1,flags=re.S)
#                selection_index_lst["Codon"] = match_content.group(1)
#                selection_index_lst["Prob[dN/dS>1]"] = match_content.group(2)
#                selection_index_lst["EBF[dN/dS]>1"] = match_content.group(3)
#                selection_index_lst["PSRF"] = match_content.group(4)
#                selection_index_lst["N_eff"] = match_content.group(5)
#                selection_index_lst["Gene_ID"] = geneid1
                match_content = re.search(r"\[RESULTS\].*\n{2}(.*)\n",records1,flags=re.S)
                if match_content is not None:
                    df = pd.read_csv(StringIO(match_content.group(1)),sep="\t")
                    df.loc[:,"Gene_ID"] = Series(geneid1,index=df.index)
                    selection_gene_df_lst1.append(df)
        if len(selection_gene_df_lst1):
            selection_dfs1 = pd.concat(selection_gene_df_lst1,axis=0,ignore_index=True)
        else:
            selection_dfs1 = pd.DataFrame(columns=['Gene_ID','values'])
        return {"screen":selection_dfs1,"csv":selection_dfs}                
    else:
        return {"csv":selection_dfs}



    
    
    
    
def parsing_MEME(p_threshold = 0.1,q_threshold = 0.05,ensembl_id = True,csv_and_screen = True):
    selection_gene_df_lst = []
    for file in glob.glob("*MEME.csv"):
        if ensembl_id == True:
            geneid = re.match(r'\S*(ENS\w{12})\S*',file).group(1)
        else:
            geneid = file
        with open(file) as f:
            df = pd.read_csv(f,sep=',',header='infer')
#            df = df.assign(Codon = Series(map(lambda x:int(x)+1,df.index.tolist()), index=df.index))
            df = df.assign(Gene_ID = Series(geneid, index=df.index))
            df = df.reset_index()
            df.rename(columns={'index':'Codon'},inplace=True)
            df['Codon'] = df['Codon']+1
            active_df = df[(df['p-value']<p_threshold) & (df['q-value']<q_threshold)]
            if len(active_df):
                selection_gene_df_lst.append(active_df)
    if len(selection_gene_df_lst):
        selection_df = pd.concat(selection_gene_df_lst,axis=0,ignore_index=True)
    else:
        selection_df = pd.DataFrame(columns=['Gene_ID','values'])
    if csv_and_screen == True:
        selection_gene_df_lst1 = []
        for file1 in glob.glob("*MEME*screen*"):
            if ensembl_id == True:
                geneid1 = re.match(r'\S*(ENS\w{12})\S*',file1).group(1)
            else:
                geneid1 = file1
            with open(file1) as f1:
                records1 = f1.read()
                match_content = re.findall(r"\|\s+(.+)\s+\*P",records1)
                if len(match_content):
                    selection_codon_dict_lst=[]
                    for line in match_content:
                        selection_codon_dict={}
                        para_list = line.split("|")
                        for para in para_list:
                            para1 = para.split(":")
                            selection_codon_dict[para1[0].strip()] = para1[1].strip()
                        selection_codon_dict_lst.append(selection_codon_dict)
                    selection_gene_df = DataFrame(selection_codon_dict_lst)
                    selection_gene_df = selection_gene_df.assign(Gene_ID = Series(geneid1, index=selection_gene_df.index))
                    selection_gene_df_lst1.append(selection_gene_df)
        if len(selection_gene_df_lst1):
            selection_df1 = pd.concat(selection_gene_df_lst1,axis=0,ignore_index=True)
            selection_df1 = selection_df1[selection_df1["p"].astype(float) < p_threshold]
        else:
            selection_df1 = pd.DataFrame(columns=['Gene_ID','values'])
        return {"screen":selection_df1,"csv":selection_df}
    else:
        return {"csv":selection_df} 

                



    
    
    
def parsing_FEL(p_threshold = 0.1,ensembl_id = True,csv_and_screen = True):
    selection_df_lst = []
    for file in glob.glob("*FEL.csv"):
        if ensembl_id == True:
            geneid = re.match(r'\S*(ENS\w{12})\S*',file).group(1)
        else:
            geneid = file
        with open(file) as f:
            df = pd.read_csv(f,sep=",",header="infer")
            df = df.reset_index()
            df.rename(columns={'index':'Codon'},inplace=True)
            df["Codon"] = df["Codon"]+1
            df = df.assign(Gene_ID = Series(geneid,index=df.index))
            active_df = df[(df["dN"] > df["dS"]) & (df["dN/dS"] > 1) & (df["p-value"] < p_threshold)]
            if len(active_df):
                selection_df_lst.append(active_df)
    if len(selection_df_lst):
        selection_df = pd.concat(selection_df_lst,axis=0,ignore_index=True)
    else:
        selection_df = pd.DataFrame(columns=['Gene_ID','values'])
    if csv_and_screen == True:
        selection_df_lst1 = []
        for file1 in glob.glob("*FEL*screen*"):
            if ensembl_id == True:
                geneid1 = re.match(r'\S*(ENS\w{12})\S*',file1).group(1)
            else:
                geneid1 = file1
            with open(file1) as f1:
                records1 = f1.read()
                match_content = re.findall(r"\|\s+(.+)\s+\*P",records1)
                if len(match_content):
                    para_dict_lst=[]
                    for line in match_content:
                        para_dict = {}
                        paras = line.split("|")
                        for para in paras:
                            para_lst = para.split(":")
                            para_dict[para_lst[0].strip()] = para_lst[1].strip()
                        para_dict_lst.append(para_dict)
                    df1 = DataFrame(para_dict_lst)
                    df1 = df1.assign(Gene_ID=Series(geneid1,index=df1.index))
                    selection_df_lst1.append(df1)
        if len(selection_df_lst1):
            selection_df1 = pd.concat(selection_df_lst1,axis=0,ignore_index=True)
            selection_df1 = selection_df1[(selection_df1["dN"].astype(float) > selection_df1["dS"].astype(float)) & (selection_df1["dN/dS"].astype(float) >1) & (selection_df1["p"].astype(float) < p_threshold)]
        else:
            selection_df1 = pd.DataFrame(columns=['Gene_ID','values'])
        return {"screen":selection_df1,"csv":selection_df}
    else:
        return {"csv":selection_df} 
            

    
        

                

def parsing_SLAC(p_threshold = 0.1,ensembl_id = True,csv_and_screen = True):
    selection_df_lst = []
    for file in glob.glob("*SLAC.csv"):
        if ensembl_id == True:
            geneid = re.match(r'\S*(ENS\w{12})\S*',file).group(1)
        else:
            geneid = file
        with open(file) as f:
            df = pd.read_csv(f,header="infer",sep='\t')
            df = df.reset_index()
            df.rename(columns={"index":"Codon"},inplace=True)
            df = df.assign(Gene_ID = Series(geneid,index=df.index))
            df["Codon"] = df["Codon"]+1
            active_df = df[(df["dN"] > df["dS"]) & (df["dN-dS"] > 0) & (df["P{S leq. observed}"] < p_threshold)]
            if len(active_df):
                selection_df_lst.append(active_df)
    if len(selection_df_lst):
        selection_df = pd.concat(selection_df_lst,axis=0,ignore_index=True)
    else:
        selection_df = pd.DataFrame(columns=['Gene_ID','values'])
    if csv_and_screen == True:
        selection_df_lst1 = []
        for file1 in glob.glob("*SLAC*screen*"):
            if ensembl_id == True:
                geneid1 = re.match(r'\S*(ENS\w{12})\S*',file1).group(1)
            else:
                geneid1 = file1
            with open(file1) as f1:
                records1 = f1.read()
                match_content = re.search(r"POSITIVELY\sSELECTED\sSITES\s\*{8}\n{3}(.+)\n{2}\*{7}\sFOUND",records1,flags=re.S)
                if match_content is not None:
#                    match_content1 = match_content.group(1).replace('-','').replace('+','').replace(' ','').replace('|',',')
                    line_lst = []
                    for idx, line in enumerate(StringIO(match_content.group(1))):
                        if idx % 2 != 0:
                            line_lst.append(line)
#                    line_lst = [re.sub('^\|', '',x) for x in line_lst]
#                    line_lst = [re.sub('\|\s*\n', '\n',x) for x in line_lst]
                    match_content1 = ''.join(line_lst).replace(' ','').replace('|',',')
                    df1 = pd.read_csv(StringIO(match_content1),sep=",")
                    df1 = df1.drop(df1.columns[[0, 5]], axis=1)
                    df1 = df1.assign(Gene_ID = Series(geneid1,index=df1.index))
                    df1 = df1[df1['dN-dS'] > 1 & (df1['p-value'] < p_threshold)]
                    if len(df1):
                        selection_df_lst1.append(df1)
        if len(selection_df_lst1):
            selection_df1 = pd.concat(selection_df_lst1,axis=0,ignore_index=True)
        else:
            selection_df1 = pd.DataFrame(columns=['Gene_ID','values'])
        return {"screen":selection_df1,"csv":selection_df}
    else:
        return {"csv":selection_df} 
        
                
                                        
if method == "s":
        df = parsing_SLAC().get("csv")
        if not df.empty:       
            df.set_index("Gene_ID").sort_index().to_csv('SLAC_parsing_json_results.csv',header=True,index=True)
        else:
            print("No codons were found to be under positive selection at p <= 0.1")
        df1 = parsing_SLAC().get("screen")
        if len(df1.index):
            df1.set_index("Gene_ID").sort_index().to_csv('SLAC_parsing_scren_results.csv',header=True,index=True)
        else:
            print("No codons were found to be under positive selection at p <= 0.1")
            
                                    
                                                
                                                                
                                                                                
                                                                                                
                                                                                                                                
            
            
if method == "fel":
        df = parsing_FEL().get("csv")
        if not df.empty:       
            df.set_index("Gene_ID").sort_index().to_csv('FEL_parsing_json_results.csv',header=True,index=True)
        else:
            print("No codons were found to be under positive selection at p <= 0.1")
        df1 = parsing_FEL().get("screen")
        if len(df1.index):
            df1.set_index("Gene_ID").sort_index().to_csv('FEL_parsing_scren_results.csv',header=True,index=True)
        else:
            print("No codons were found to be under positive selection at p <= 0.1")
    
    


                                                                
                                                                
                                                                                                
                                                                                                                                                                
                                                                
if method == "m":
    with pd.option_context('display.precision',5):
        df = parsing_MEME().get("csv")
        if not df.empty:       
            df.set_index("Gene_ID").sort_index().to_csv('MEME_parsing_json_results.csv',header=True,index=True)
        else:
            print("No codons were found to be under positive selection at p <= 0.1")
        df1 = parsing_MEME().get("screen")
        if len(df1.index):
            df1.set_index("Gene_ID").sort_index().to_csv('MEME_parsing_scren_results.csv',header=True,index=True)
        else:
            print("No codons were found to be under positive selection at p <= 0.1")
                    
                
        
    


if method == "f":
    with pd.option_context('display.precision',5):
        df = parsing_FUBAR().get("csv")
        if not df.empty:       
            df.set_index("Gene_ID").sort_index().to_csv('FUBAR_parsing_json_results.csv',header=True,index=True)
        else:
            print("No codons were found to be under positive selection at prob >= 0.9")
        df1 = parsing_FUBAR().get("screen")
        if len(df1.index):
            df1.set_index("Gene_ID").sort_index().to_csv('FUBAR_parsing_scren_results.csv',header=True,index=True)
        else:
            print("No codons were found to be under positive selection at prob >= 0.9")
    

                            
                
            
                                                    
                                                            
                                                                                    
        
if method == "r":
    with pd.option_context('display.precision',5):
#    pd.options.display.float_format = '{:,.5f}'.format
        df = parsing_RELAX().get("json")
        if not df.empty:       
            df.set_index("Gene_ID").sort_index().to_csv('RELAX_parsing_json_results.csv',header=True,index=True)
        else:
            print("No branches were found to be under relaxed selection at p <= 0.05")
        df1 = parsing_RELAX().get("screen")
        if len(df1.index):
            df1.set_index("Gene_ID").sort_index().to_csv('RELAX_parsing_scren_results.csv',header=True,index=True)
        else:
            print("No branches were found to be under relaxed selection at p <= 0.05")
               
            
        





if method == 'b':
    with pd.option_context('display.precision',5):
        df = parsing_BUSTED().get("json_foreground")
        if len(df):
            df.set_index("Gene_ID").sort_index().to_csv('BUSTED_parsing_json_foreground.csv',header=True,index=True)
        else:
            print("No branches were found to be under positive selection at p <= 0.05")
        df2 = parsing_BUSTED().get("json_full")
        if len(df2):
            df2.set_index("Gene_ID").sort_index().to_csv('BUSTED_parsing_json_full.csv',header=True,index=True)
        else:
            print("No json_full or No branches were found to be under positive selection at p <= 0.05")
        df1 = parsing_BUSTED().get("screen")
        if len(df1):
            df1.set_index("Gene_ID").sort_index().to_csv('BUSTED_parsing_scren_results.csv',header=True,index=True)
        else:
            print("No branches were found to be under positive selection at p <= 0.05")
           
            


"""   
if method == 'a':    
    with open(method_dict[method]+"_parsing_results.jsons", "w") as op:
        json.dump(parsing_aBSREL()["json"],op)
    with open(method_dict[method]+"_parsing_results_scren.txt", "w") as op:
        op.write(json.dumps(parsing_aBSREL()["screen"]))
    with open(method_dict[method]+"_parsing_results.json.txt", "w") as op:
        for gene in parsing_aBSREL()["json"]:
            op.write(gene+'\n')
            for sp in parsing_aBSREL()["json"][gene]:
                op.write(sp+'\t'+str(parsing_aBSREL()["json"][gene][sp]["LRT"])+'\t'+str(parsing_aBSREL()["json"][gene][sp]["p"])+'\n')
    with open(method_dict[method]+"_parsing_results.json.geneid.txt", "w") as op:
        for gene in parsing_aBSREL()["json"]:
            op.write(gene+'\n')
    with open("aBSREL_parsing_results.json1.txt", "w") as op:
        for gene in parsing_aBSREL()["json"]:
#            op.write(gene+'\t')
            for sp in parsing_aBSREL()["json"][gene]:
                op.write(gene+'\t'+sp+'\t'+str(parsing_aBSREL()["json"][gene][sp]["LRT"])+'\t'+str(parsing_aBSREL()["json"][gene][sp]["p"])+'\n')
    aBSREL_table = pd.read_table("aBSREL_parsing_results.json1.txt",sep="\t",header=None,names=["Gene_ID","Taxon","LRT","p-value"])
    aBSREL_table.to_csv("aBSREL_parsing_results.json1.csv",header=True,index=False)
"""
if method == 'a':
    with pd.option_context('display.precision',5):
        df = parsing_aBSREL().get("json")
        if not df.empty:       
            df.set_index("Gene_ID").sort_index().to_csv('aBSREL_parsing_json_results.csv',header=True,index=True)
        else:
            print("No branches were found to be under positive selection at p <= 0.05")
        df1 = parsing_aBSREL().get("screen")
        if len(df1.index):
            df1.set_index("Gene_ID").sort_index().to_csv('aBSREL_parsing_scren_results.csv',header=True,index=True)
        else:
            print("No branches were found to be under positive selection at p <= 0.05")


       

