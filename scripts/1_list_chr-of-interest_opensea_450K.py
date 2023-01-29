#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


# ## ToDo: Specify target chromosome and subregion
# - chromosome: 1 ~ 22
# - subregion: island / opensea / shore / shelf

# ## Data
# - 450K: GSE36369 (EBV)
# - Hi-C: GSE18199 (HiC-EBV-2009)

# # Header
# - Chromosome_36: Chromosome genome build 36
# - Coordinate_36: Coordinates genome build 36
# - Relation_to_UCSC_CpG_Island: Relationship to Canonical CpG Island: Shores - 0-2 kb from CpG island; Shelves - 2-4 kb from CpG island.
#     - ['N_Shore' 'S_Shelf' nan 'Island' 'S_Shore' 'N_Shelf']

# In[4]:

suppl_fname = 'GPL13534_HumanMethylation450_15017482_v.1.1.csv'
manifest_fname = 'GPL13534_450K_Manifest_header_Descriptions.xlsx'

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', required=True, default = '14', help = 'chromosome-of-interest', type = str)
    parser.add_argument('-w_dir', '--working_dir', required = True, help = 'working directory', default = '/data/project/jeewon/research/3dith-reproduce', type = str)
    parser.add_argument('-d_dir', '--data_dir', required = True, help = '450K data directory', default = '/data/project/jeewon/research/3dith-reproduce/data', type = str)
    parser.add_argument('-s_dir', '--save_dir', required = True, help = 'saving directory', default = '/data/project/jeewon/research/3dith-reproduce/result', type = str)
    parser.add_argument('--cpg_type', required = True, help = 'island, opensea, shelf, shore', type = str)
    parser.add_argument('--grch', required = True, help = 'GRCh version', default = 36, type = int)
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    if not os.path.exists(args.data_dir):
        os.makedirs(args.data_dir)
    if not os.path.exists(args.working_dir):
        os.makedirs(args.working_dir)
    
    os.chdir(args.working_dir)
    
    if suppl_fname in os.listdir(args.data_dir):
        suppl = pd.read_csv(os.path.join(args.data_dir, suppl_fname), skiprows = [0, 1, 2, 3, 4, 5, 6], index_col = 'IlmnID')
    else:
        os.chdir(args.data_dir)
        cmd1 = 'wget https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/suppl/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz'
        cmd2 = 'gzip -d GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz'
        os.system(cmd1)
        os.system(cmd2)
        os.chdir(args.working_dir)
        suppl = pd.read_csv(os.path.join(args.data_dir, suppl_fname), skiprows = [0, 1, 2, 3, 4, 5, 6], index_col = 'IlmnID')
    
    '''
    if  manifest_fname in os.listdir(args.data_dir):
        manifest = pd.read_excel(os.path.join(args.data_dir, manifest_fname))
    else:
        os.chdir(args.data_dir)
        cmd1 = 'wget https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/suppl/GPL13534_450K_Manifest_header_Descriptions.xlsx.gz''
        cmd2 = 'gzip -d GPL13534_450K_Manifest_header_Descriptions.xlsx.gz'
        os.system(cmd1)
        os.system(cmd2)
        os.chdir(args.working_dir)
        manifest = pd.read_excel(os.path.join(args.data_dir, manifest_fname))
    '''
    # extract hg18 info from suppl.
    suppl2 = suppl[['Chromosome_36', 'Coordinate_36', 'Relation_to_UCSC_CpG_Island']].copy()
    
    # extract cg probes
    suppl3 = suppl2.loc[['cg' in suppl2.index.values[i] for i in range(len(suppl2.index.values))],:].copy()
    
    # extract chromosome-of-interest
    suppl4 = suppl3[suppl3['Chromosome_36'].isin([args.chrom, int(args.chrom)])].copy()
    
    if args.cpg_type=='opensea':
        suppl5 = suppl4[suppl4['Relation_to_UCSC_CpG_Island'].isna()].copy()
    elif args.cpg_type=='shelf':
        shelf_mask = ['Shelf' in str(suppl4['Relation_to_UCSC_CpG_Island'].values[i]) for i in range(suppl4.shape[0])]
        suppl5 = suppl4.iloc[shelf_mask,:].copy()
    elif args.cpg_type=='shore':
        shore_mask = ['Shore' in str(suppl4['Relation_to_UCSC_CpG_Island'].values[i]) for i in range(suppl4.shape[0])]
        suppl5 = suppl4.iloc[shore_mask,:].copy()
    elif args.cpg_type=='island':
        island_mask = ['Island' in str(suppl4['Relation_to_UCSC_CpG_Island'].values[i]) for i in range(suppl4.shape[0])]
        suppl5 = suppl4.iloc[island_mask,:].copy()
    else:
        print("Invalid region! Region should be one of ['opensea', 'shelf', 'shore', 'island'].")
        
    # save result
    suppl5_fname = f'450K-GRCh{args.grch}-chr{args.chrom}-{args.cpg_type}.csv'
    suppl5.to_csv(os.path.join(args.save_dir, suppl5_fname))
    print(os.path.join(args.save_dir, suppl5_fname))

