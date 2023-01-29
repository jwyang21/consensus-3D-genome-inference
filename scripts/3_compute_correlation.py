#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import argparse


# In[2]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', required=True, default = '14', help = 'chromosome-of-interest', type = str)
    parser.add_argument('-w_dir', '--working_dir', required = True, help = 'working directory', default = '/data/project/jeewon/research/3dith-reproduce', type = str)
    parser.add_argument('-d_dir', '--data_dir', required = True, help = '450K data directory', default = '/data/project/jeewon/research/3dith-reproduce/data', type = str)
    parser.add_argument('-s_dir', '--save_dir', required = True, help = 'saving directory', default = '/data/project/jeewon/research/3dith-reproduce/result', type = str)
    parser.add_argument('--cpg_type', required = True, help = 'island, opensea, shelf, shore', type = str)
    parser.add_argument('--grch', required = True, help = 'GRCh version', default = 36, type = int)
    return parser.parse_args()


# In[3]:


suppl_fname = 'GPL13534_HumanMethylation450_15017482_v.1.1.csv'
manifest_fname = 'GPL13534_450K_Manifest_header_Descriptions.xlsx'
beta_fname = 'GSE36369_series_matrix.txt'


# In[ ]:


if __name__ == '__main__':

    args = parse_arguments()
    if not os.path.exists(args.data_dir):
        os.makedirs(args.data_dir)
    if not os.path.exists(args.working_dir):
        os.makedirs(args.working_dir) 
    
    os.chdir(args.working_dir)
    
    # 1. Load beta and cpg probe data
    beta_fname = os.path.join(args.save_dir, f'450K-GRCh{args.grch}-chr{args.chrom}-{args.cpg_type}-beta.csv')
    cpg_fname = os.path.join(args.save_dir, f'450K-GRCh{args.grch}-chr{args.chrom}-{args.cpg_type}.csv')
    beta = pd.read_csv(beta_fname, index_col = 0)
    cpg = pd.read_csv(cpg_fname, index_col = 0)
    
    # 2. Delete columns with more than 100 NaNs.
    beta2 = beta.iloc[:,beta.isna().sum().values <= 100].copy()
    beta3 = beta2.iloc[(beta2.T.isna().sum().values == 0),:].copy() # (row: CpG probe, column: sample)
    beta3_fname = os.path.join(args.save_dir, f'450K-GRCh{args.grch}-chr{args.chrom}-{args.cpg_type}-beta.cleaned.csv')
    beta3.to_csv(beta3_fname)
    print(beta3_fname)
    
    # 3. Calculate correlation matrix.
    corr = np.zeros((beta3.shape[0], beta3.shape[0]), dtype = float)
    for i in range(beta3.shape[0]): #각각의 CpG probe에 대해
        if i % (beta3.shape[0]//10) == 0:
            print("\n===")
            print("Processing {}-th cg probe.".format(i))
        current_cg = beta3.iloc[i,:].values.flatten()
        for j in range(beta3.shape[0]): #현재 보고 있는 CpG probe - 모든 CpG probe 간 all-pairwise PCC 계산
            current_probe_cg = beta3.iloc[j,:].values.flatten()
            pearson = stats.pearsonr(current_cg, current_probe_cg)[0]
            corr[i,j] = pearson
    corr_df = pd.DataFrame(corr, index = beta3.index.values, columns = beta3.index.values)
    corr_df_fname = os.path.join(args.save_dir, f'450K-GRCh{args.grch}-chr{args.chrom}-{args.cpg_type}-corr.csv')
    corr_df.to_csv(corr_df_fname)
    print(corr_df_fname)
    assert corr_df.shape[0] == beta3.shape[0] #number of total CpG probes

