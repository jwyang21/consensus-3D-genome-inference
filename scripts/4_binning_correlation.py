#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import argparse


# In[ ]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', required=True, default = '14', help = 'chromosome-of-interest', type = str)
    parser.add_argument('-w_dir', '--working_dir', required = True, help = 'working directory', default = '/data/project/jeewon/research/3dith-reproduce', type = str)
    parser.add_argument('-d_dir', '--data_dir', required = True, help = '450K data directory', default = '/data/project/jeewon/research/3dith-reproduce/data', type = str)
    parser.add_argument('-s_dir', '--save_dir', required = True, help = 'saving directory', default = '/data/project/jeewon/research/3dith-reproduce/result', type = str)
    parser.add_argument('--cpg_type', required = True, help = 'island, opensea, shelf, shore', type = str)
    parser.add_argument('--grch', required = True, help = 'GRCh version', default = 36, type = int)
    parser.add_argument('--chrom_len', required = True, help = 'chromosome length filename', default = '/data/project/jeewon/research/reference/hg18.fa.sizes', type = str)
    parser.add_argument('--binsize', required = True, help = 'binsize of binned correlation matrix', default = int(1e5), type = int)
    return parser.parse_args()


# ### ToDo
# #### 1. Make CpG dictionary
# - key: bin name (chr14:[start_index]-[end_index]
# - dictionary: name of cg probes in this bin
# - Find out which open sea CpG probes are in each genomic bin of fixed binsize.   
# 
# #### 2. Binning the raw corr. matrix
# - Only consider bins containing open sea CpG probes
# - ex) All-pairwise Pearson correlation values between {beta value of CpG probes in bin 1} and {beta value of CpG probes in bin 2} -> compute median -> correlation between bin 1 and bin 2

# In[ ]:


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
    
    # 1. Load data
    ## 1-1. Import 450K-derived correlation matrix
    corr_df_fname = os.path.join(args.save_dir, f'450K-GRCh{args.grch}-chr{args.chrom}-{args.cpg_type}-corr.csv')
    corr = pd.read_csv(corr_df_fname, index_col = 0)
    
    ## 1-2. CpG list
    cpg_fname = os.path.join(args.save_dir, f'450K-GRCh{args.grch}-chr{args.chrom}-{args.cpg_type}.csv')
    cpg = pd.read_csv(cpg_fname, index_col = 0)
    cpg2 = cpg[cpg.index.isin(corr.index.values)].copy() # Extract cg probes included in corr. matrix
    
    ## 1-3. chromosome length (to bin the correlation matrix)
    chr_length = pd.read_csv(args.chrom_len, sep = '\t', header = None, names = ['chr','length'])
    chr_length = chr_length.astype({'length':'int'})
    chr_interest_length = chr_length[chr_length['chr']==str('chr'+args.chrom)].length.values[0]
    print(f"chr{args.chrom} length: {chr_interest_length}")
    
    # 2. Compute binned correlation matrix
    ## 2-1. Bin the cpg list
    bin_num = (chr_interest_length // args.binsize)+1
    print(f"binsize: {args.binsize}, number of bins: {bin_num}")
    
    bin_probe = {}
    bin_with_probes = []
    bin_with_probes_index = []
    for i in range(bin_num):
        if i % (bin_num//10) == 0:
            print("Processing {}-th bin".format(i))
        bin_name = 'chr'+args.chrom+':'+str(args.binsize*i)+'-'+str(args.binsize*(i+1)-1) #[start, end]
        bin_start = args.binsize * i
        bin_end = args.binsize * (i+1) - 1
        tmp = [] # list of CpG probes located in this bin
        for j in range(cpg2.shape[0]): 
            current_cg = cpg2.index.values[j]
            current_cg_pos = int(cpg2['Coordinate_36'].values[j])
            if current_cg_pos >= bin_start and current_cg_pos <= bin_end:
                tmp.append(current_cg)
        bin_probe[bin_name] = tmp
        if len(tmp)>0:
            bin_with_probes.append(bin_name) 
            bin_with_probes_index.append(i) 
            
    ## 2-2. save binned CpG probe
    npz_fname = os.path.join(args.save_dir, f'bin_probe_GRCh{args.grch}_chr{args.chrom}_{args.cpg_type}')
    np.savez(npz_fname, **bin_probe)
    
    bin_with_probes = np.array(bin_with_probes)
    bin_with_probes_npy = os.path.join(args.save_dir, f'bin_with_probes_GRCh{args.grch}_chr{args.chrom}_{args.cpg_type}')
    np.save(bin_with_probes_npy, bin_with_probes)
    
    bin_with_probes_index = np.array(bin_with_probes_index)
    bin_with_probes_index_npy = os.path.join(args.save_dir, f'bin_with_probes_index_GRCh{args.grch}_chr{args.chrom}_{args.cpg_type}')
    np.save(bin_with_probes_index_npy, bin_with_probes_index)

    ## 2-3. Bin the correlation matrix
    print("\nBinning the raw corr. matrix")
    binned_corr = np.zeros((len(bin_with_probes),len(bin_with_probes)), dtype = float)
    for i in range(len(bin_with_probes)):
        if i % (len(bin_with_probes)//10) == 0:
            print("\n===============================")
            print("Processing {}-th bin_with_probe".format(i))
        current_bin = bin_with_probes[i]
        current_bin_cg = bin_probe[current_bin]
        for j in range(len(bin_with_probes)):
            current_bin2 = bin_with_probes[j]
            current_bin2_cg = bin_probe[current_bin2]
            tmp = []
            for k in range(len(current_bin_cg)):
                for l in range(len(current_bin2_cg)):
                    tmp.append(corr.loc[current_bin_cg[k]][current_bin2_cg[l]])
            binned_corr[i,j] = np.median(tmp)
    binned_corr_df = pd.DataFrame(binned_corr, index = bin_with_probes, columns = bin_with_probes)
    binned_corr_df_fname = os.path.join(args.save_dir, f'450K-GRCh{args.grch}-chr{args.chrom}-{args.cpg_type}-corr.binned.csv')
    binned_corr_df.to_csv(binned_corr_df_fname)
    print(binned_corr_df_fname)

