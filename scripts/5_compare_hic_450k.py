#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
import argparse


# ## README
# - input: binned correlation matrices made from Hi-C and 450K, respectively.
#     - 450K: GSE36369 (EBV)
#     - Hi-C: GSE18199 (HiC-EBV-2009)
# - output: 
#     - Heatmap of each binned correlation matrix.
#     - A/B compartment estimated from each binned correlation matrix.

# ### PCA
# - to compute A/B compartment distribution
# - binned_corr_hic2 // binned_corr_450k
# - Extract PC1
# - get sign of each entry of PC1
# - assign each sign (-/+) to A/B compartment
# - calculate aggrement ratio of A/B compartment calculated from HiC and 450K data.
# - Exclude bins which satisfies abs(entry of PC1 of this bin) <= 0.01 --> thresholding

# ### Assign A/B compartments to each bin in binned corr. matrices
# - The sign of the eigenvector is chosen so that the sign of the correlation between the eigenvector and column sums of the correlation matrix is positive; this ensures that positive values of the eigenvector are associated with the closed compartment
# - binned_corr_450k
# - binned_corr_hic2
# - pca_hic_pc1
# - pca_450k_pc1

# ### Densities of corr. of the 450k methylation probes
# - binned_corr_450k
# - binned_corr_hic2
# - ===================
# - pca_hic_pc1
# - pca_450k_pc1_v2
# - ===================
# - pca_hic_thresholded
# - pca_450k_v2_thresholded

# #### positive values of the eigenvector are associated with the closed compartment
# - (+) --> B (closed chromatin) // (-) --> A (open chromatin)

# In[1]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', required=True, default = '14', help = 'chromosome-of-interest', type = str)
    parser.add_argument('-w_dir', '--working_dir', required = True, help = 'working directory', default = '/data/project/jeewon/research/3dith-reproduce', type = str)
    parser.add_argument('-d_dir', '--data_dir', required = True, help = '450K data directory', default = '/data/project/jeewon/research/3dith-reproduce/data', type = str)
    parser.add_argument('-s_dir', '--save_dir', required = True, help = 'saving directory', default = '/data/project/jeewon/research/3dith-reproduce/result', type = str)
    parser.add_argument('--cpg_type', required = True, help = 'island, opensea, shelf, shore', type = str)
    parser.add_argument('--grch', required = True, help = 'GRCh version', default = 36, type = int)
    return parser.parse_args()


# In[ ]:


suppl_fname = 'GPL13534_HumanMethylation450_15017482_v.1.1.csv'
manifest_fname = 'GPL13534_450K_Manifest_header_Descriptions.xlsx'
beta_fname = 'GSE36369_series_matrix.txt'


# In[1]:


'''
chrom, working_dir, data_dir, save_dir, cpg_type, grch = 14, '/data/project/jeewon/research/3dith-reproduce', '/data/project/jeewon/research/3dith-reproduce/data', \
'/data/project/jeewon/research/3dith-reproduce/result', 'opensea', 36
'''


# In[ ]:


if __name__ == '__main__':

    args = parse_arguments()
    if not os.path.exists(args.data_dir):
        os.makedirs(args.data_dir)
    if not os.path.exists(args.working_dir):
        os.makedirs(args.working_dir) 
    
    os.chdir(args.working_dir)
    
    # 1. Import correlation matrix
    ## 1-1. Import 450K binned correlation matrix
    binned_corr_450k_fname = os.path.join(args.save_dir, f'450K-GRCh{args.grch}-chr{args.chrom}-{args.cpg_type}-corr.binned.csv')
    binned_corr_450k = pd.read_csv(binned_corr_450k_fname, index_col = 0)
    
    ## 1-2. Import HiC binned correlation matrix
    binned_corr_hic_fname = os.path.join(args.data_dir, f'HIC_gm06690_chr{args.chrom}_chr{args.chrom}_100000_pearson.txt')
    if not os.path.exists(binned_corr_hic_fname):
        os.chdir(args.data_dir)
        cmd1 = 'wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18199/suppl/GSE18199_binned_heatmaps.zip.gz'
        cmd2 = 'gzip -d GSE18199_binned_heatmaps.zip.gz'
        cmd3 = 'unzip GSE18199_binned_heatmaps.zip'
        cmd4 = 'rm -rf HIC_*exp.txt'
        cmd5 = 'rm -rf HIC_*obs.txt'
        os.system(cmd1)
        os.system(cmd2)
        os.system(cmd3)
        os.system(cmd4)
        os.system(cmd5)
        os.chdir(args.working_dir)
    binned_corr_hic = pd.read_csv(os.path.join(args.data_dir, binned_corr_hic_fname), skiprows=[0], index_col = 0, sep = '\t')
    binned_corr_hic = binned_corr_hic.iloc[:,:-1].copy()
    
    ## 1-3. Import CpG information
    
    #bin_with_probes_list_npy = os.path.join(args.save_dir, f'bin_with_probes_GRCh{args.grch}_chr{args.chrom}_{args.cpg_type}.npy')
    #bin_with_probes_list = np.load(bin_with_probes_list_npy)
    
    bin_with_probes_index_npy = os.path.join(args.save_dir, f'bin_with_probes_index_GRCh{args.grch}_chr{args.chrom}_{args.cpg_type}.npy')
    bin_with_probes_index = np.load(bin_with_probes_index_npy)
    
    ## 1-4. Extract target bins from binned HiC correlation matrix
    binned_corr_hic2 = binned_corr_hic.iloc[bin_with_probes_index, bin_with_probes_index].copy()
    
    # 2. Plot (check whether plots are consistent in Hi-C and 450K)
    ## 2-1. Histogram of binned correlation values from Hi-C and 450K
    fig = plt.figure(figsize = (7,3.5))
    ax1 = fig.add_subplot(121)
    ax1.hist(binned_corr_450k.values.flatten())
    ax1.set_title('Histogram of binned corrmat (450K)')
    ax2 = fig.add_subplot(122)
    ax2.hist(binned_corr_hic2.values.flatten())
    ax2.set_title('Histogram of binned corrmat (Hi-C)')
    suptitle_ = f'Histograms of binned corrmats (chr{args.chrom}, {args.cpg_type})'
    plt.suptitle(suptitle_, fontsize = 15)
    plt.tight_layout()
    #plt.show()
    hist_fname = os.path.join(args.save_dir, f'Hist_binned_corrmat_GRCh{args.grch}_chr{args.chrom}_{args.cpg_type}.png')
    fig.tight_layout()
    plt.savefig(hist_fname)
    print("hist_fname: ", hist_fname)
    print("\n")
    plt.clf()
    
    ## 2-2. Heatmap of binned corrmats from 450K and HiC
    fig = plt.figure(figsize = (7,3.5))
    ax1 = fig.add_subplot(121)
    #ax1.matshow(binned_corr_450k.values, cmap = 'Reds', vmin = np.min(binned_corr_450k.values.flatten()), vmax = np.max(binned_corr_450k.values.flatten()))
    ax1.matshow(binned_corr_450k.values, cmap = 'Reds')
    ax1.set_title('Heatmap of binned corrmat (450K)')
    ax2 = fig.add_subplot(122)
    #ax2.matshow(binned_corr_hic2.values, cmap = 'Reds', vmin = np.min(binned_corr_hic2.values.flatten()), vmax = np.max(binned_corr_hic2.values.flatten()))
    ax2.matshow(binned_corr_hic2.values, cmap = 'Reds')
    ax2.set_title('Heatmap of binned corrmat (Hi-C)')
    suptitle_ = f'Histograms of binned corrmats (chr{args.chrom}, {args.cpg_type})'
    plt.suptitle(suptitle_, fontsize = 15)
    plt.tight_layout()
    #plt.show()
    heatmap_fname = os.path.join(args.save_dir, f'Heatmap_binned_corrmat_GRCh{args.grch}_chr{args.chrom}_{args.cpg_type}.png')
    fig.tight_layout()
    plt.savefig(heatmap_fname)
    print("heatmap_fname: ", heatmap_fname)
    print("\n")
    plt.clf()
    
    # 3. Compare Hi-C PC1 and 450K PC1
    ## 3-1. PCA
    pca = PCA(n_components = 3)
    pca_hic = pca.fit_transform(binned_corr_hic2.values)
    
    pca2 = PCA(n_components = 3)
    pca_450k = pca2.fit_transform(binned_corr_450k.values)
    
    assert pca_hic.shape[0] == pca_450k.shape[0]
    
    pca_450k_pc1 = pca_450k[:,0]
    pca_hic_pc1 = pca_hic[:,0]
    
    ## 3-2. Assign A/B compartments to each bin in binned corr. matrices
    ### 3-2-1. Make sure that column sums of binned correlation matrix and PC1 from 450K binned correlation matrix have positive correlation.
    corr_colsum = binned_corr_450k.sum().values
    sgn = np.sign(stats.pearsonr(pca_450k_pc1, corr_colsum)[0])
    
    if sgn == -1:
        print("Multiply (-1) to 450k PC1 vector.")
        pca_450k_pc1_v2 = (-1) * pca_450k_pc1
    elif sgn == 1:
        pca_450k_pc1_v2 = pca_450k_pc1
    else: #sign == 0
        print("PCC(PC1, column sum) of 450K binned corrmat is 0. Error might have occurred.")

    ### 3-2-2. unthresholded PCC(HiC PC1, 450k PC1)
    print("Unthresholded PCC(HiC PC1, 450k PC1):", end = ' ')
    print(stats.pearsonr(pca_hic_pc1, pca_450k_pc1_v2))
    print("\n")
    
    ### 3-2-3. Scatter(HiC PC1, 450k PC1)
    fig = plt.figure(figsize = (3.5,3.5))
    ax = fig.add_subplot(111)
    ax.scatter(pca_hic_pc1, pca_450k_pc1_v2)
    ax.set_xlabel("PC1 from Hi-C")
    if sgn == -1:
        ylabel = '(-1) * PC1 from 450K'
    elif sgn == 1:
        ylabel = 'PC1 from 450K'
    ax.set_ylabel(ylabel)
    scatter_title_ = f'Scatter of PC1 (chr{args.chrom}, {args.cpg_type})'
    ax.set_title(scatter_title_)
    #plt.show()
    scatter_fname = os.path.join(args.save_dir, f'Scatter_PC1_GRCh{args.grch}_chr{args.chrom}_{args.cpg_type}.png')
    fig.tight_layout()
    plt.savefig(scatter_fname)
    print("scatter_fname: ", scatter_fname)
    print("\n")
    plt.clf()
    
    ### 3-2-4. Plot Hi-C PC1 and 450K PC1
    fig = plt.figure(figsize = (7, 3.5))

    ax1 = fig.add_subplot(211)
    ax1.plot(pca_hic_pc1)
    ax1.set_title("PC1 from Hi-C")

    ax2 = fig.add_subplot(212)
    ax2.plot(pca_450k_pc1_v2)

    if sgn == -1:
        title = '(-1) * PC1 from 450K'
    elif sgn == 1:
        title = 'PC1 from 450K'
    ax2.set_title(title)

    suptitle_ = f'PC1 (chr{args.chrom}, {args.cpg_type})'
    plt.suptitle(suptitle_ , fontsize = 15)
    plt.tight_layout()
    #plt.show()
    plot_fname = os.path.join(args.save_dir, f'PC1_chr_{args.chrom}_{args.cpg_type}.png')
    fig.tight_layout()
    plt.savefig(plot_fname)
    print("plot_fname: ", plot_fname)
    print("\n")
    plt.clf()
    
    ### 3-2-5. Thresholding
    #### Leave bins whose absolute values >= 1e-2
    threshold_mask = np.array([abs(pca_hic_pc1[i]) >= 0.01 and abs(pca_450k_pc1_v2[i]) >= 0.01 for i in range(len(pca_hic_pc1))])
    pca_450k_v2_thresholded = pca_450k_pc1_v2[threshold_mask].copy()
    pca_hic_thresholded = pca_hic_pc1[threshold_mask].copy()
    
    ### 3-2-6. Thresholded PCC (450K PC1, Hi-C PC1)
    print("Thresholded PCC(HiC PC1, 450k PC1): ")
    print(stats.pearsonr(pca_hic_thresholded, pca_450k_v2_thresholded))
    
    
    ### 3-2-6. Scatter (thresholded Hi-C PC1, thresholded 450K PC1)
    fig = plt.figure(figsize = (3.5,3.5))
    ax = fig.add_subplot(111)
    ax.scatter(pca_hic_thresholded, pca_450k_v2_thresholded)
    ax.set_xlabel("PC1 from Hi-C")

    if sgn == -1:
        ylabel = '(-1) * PC1 from 450K'
    elif sgn == 1:
        ylabel = 'PC1 from 450K'
    ax.set_ylabel(ylabel)

    title_ = f'Scatter of PC1 (chr{args.chrom}, {args.cpg_type})'
    ax.set_title(title_)
    #plt.show()
    threshold_plot_fname = os.path.join(args.save_dir, f'Scatter_thresholded_PC1_GRCh{args.grch}_chr{args.chrom}_{args.cpg_type}.png')
    fig.tight_layout()
    plt.savefig(threshold_plot_fname)
    print(threshold_plot_fname)
    plt.clf()
    
    ### 3-2-7. Plot thresholded Hi-C PC1 and 450K PC1
    fig = plt.figure(figsize = (7, 3.5))

    ax1 = fig.add_subplot(211)
    ax1.plot(pca_hic_thresholded)
    ax1.set_title("PC1 from Hi-C (thresholded)")
    ax2 = fig.add_subplot(212)
    ax2.plot(pca_450k_v2_thresholded)

    if sgn == -1:
        title = '(-1) * PC1 from 450K (thresholded)'
    elif sgn == 1:
        title = 'PC1 from 450K (thresholded)'
    ax2.set_title(title)

    suptitle_ = f'Thresholded PC1 (chr{args.chrom}, {args.cpg_type})'
    plt.suptitle(suptitle_, fontsize = 15)
    plt.tight_layout()
    #plt.show()
    plot_fname = os.path.join(args.save_dir, f'PC1_thresholded_chr_{args.chrom}_{args.cpg_type}.png')
    fig.tight_layout()
    plt.savefig(plot_fname)
    print("threshold_plot_fname: ", plot_fname)
    print("\n")
    plt.clf()
    
    # 4. Compute A/B compartment distribution
    ## 4-1. Agreement of A/B compartment distribution
    comp_450k = ["A" if np.sign(x)==-1 else "B" for x in pca_450k_pc1_v2]
    comp_hic = ["A" if np.sign(x)==-1 else "B" for x in pca_hic_pc1]
    
    print("A/B agreement between 450K and Hi-C: ", end = ' ')
    print(np.array([comp_450k[i]==comp_hic[i] for i in range(len(comp_hic))]).sum() / len(comp_hic))
    print("\n")
    
    # 5. Density plot (open-open, open-closed, closed-closed)
    ## 5-1. Compute open-open, open-closed, closed-closed correlation in 450K and HiC binned corrmat.
    open_open_450k = []
    open_closed_450k = []
    closed_closed_450k = []

    open_open_hic = []
    open_closed_hic = []
    closed_closed_hic = []

    for i in range(len(comp_450k)): #since len(comp_450k) == len(comp_hic), doesn't matter whatever you use.
        if i % (len(comp_450k) // 10) == 0:
            print("============================================================")
            print("Processing {}-th bin".format(i))
        bin_comp_450k = comp_450k[i]
        bin_comp_hic = comp_hic[i]
        for j in range(len(comp_450k)):
            if i % (len(comp_450k) // 10) == 0 and j % (len(comp_450k) // 10) == 0:
                print("\nProcessing {}-th probe bin".format(j))

            bin2_comp_450k = comp_450k[j]
            bin2_comp_hic = comp_hic[j]
            if bin_comp_450k == bin2_comp_450k:
                if bin_comp_450k == 'A': #open - open
                    open_open_450k.append(binned_corr_450k.iloc[i][j])
                else: #bin_comp_450k == 'B' #closed = closed
                    closed_closed_450k.append(binned_corr_450k.iloc[i][j])
            if bin_comp_450k != bin2_comp_450k: #open - closed
                open_closed_450k.append(binned_corr_450k.iloc[i][j])
            if bin_comp_hic == bin2_comp_hic:
                if bin_comp_hic == 'A': #open - open
                    open_open_hic.append(binned_corr_hic2.iloc[i][j])
                else: #bin_comp_hic == 'B' #closed = closed
                    closed_closed_hic.append(binned_corr_hic2.iloc[i][j])
            if bin_comp_hic != bin2_comp_hic: #open - closed
                open_closed_hic.append(binned_corr_hic2.iloc[i][j])
                
    print("Number of entries from each compartment-compartment binned correlations: ")
    print("open_open_450k: ", end = ' ')
    print(len(open_open_450k))
    print("open_closed_450k: ", end = ' ')
    print(len(open_closed_450k))
    print("closed_closed_450k: ", end = ' ')
    print(len(closed_closed_450k))
    print("total 450k: ", end = ' ')
    print(len(binned_corr_450k.values.flatten()))
    
    assert len(open_open_450k) + len(open_closed_450k) + len(closed_closed_450k) - len(binned_corr_450k.values.flatten()) == 0
    # open_open_450k + open_closed_450k + closed_closed_450k - total 450k (should be 0):  0
    
    print("open_open_hic: ", end = ' ')
    print(len(open_open_hic))
    print("open_closed_hic: ", end = ' ')
    print(len(open_closed_hic))
    print("closed_closed_hic: ", end = ' ')
    print(len(closed_closed_hic))
    print("total hic: ", end = ' ')
    print(len(binned_corr_hic2.values.flatten()))
    print("\n")
    
    assert len(open_open_hic) + len(open_closed_hic) + len(closed_closed_hic) - len(binned_corr_hic2.values.flatten()) == 0
    # open_open_hic + open_closed_hic + closed_closed_hic - total hic(should be 0):  0
    
    ## 5-2. Plot using KDE (kernel density estimation)
    ### To capture the structure of long-range correlations in DNA methylation data
    ### Characteristic: The lack of decay of corr. with distance
    ### refer to https://gist.github.com/mwaskom/de44147ed2974457ad6372750bbe5751
    
    fig = plt.figure(figsize = (7,3.5))
    ax1 = fig.add_subplot(121)
    #sns.distplot(open_open_450k, kde=True, hist=True, ax = ax1, label = 'Open - Open (450K)')
    #sns.distplot(open_closed_450k, kde = True, hist=True, ax = ax1, label = 'Open - Closed (450K)')
    #sns.distplot(closed_closed_450k, kde = True, hist=True, ax = ax1, label = 'Closed - Closed (450K)')
    sns.histplot(open_open_450k, kde=True, ax = ax1, label = 'Open - Open (450K)')
    sns.histplot(open_closed_450k, kde = True, ax = ax1, label = 'Open - Closed (450K)')
    sns.histplot(closed_closed_450k, kde = True, ax = ax1, label = 'Closed - Closed (450K)')
    ax1.set_title("Binned correlations (450K)")# Gaussian KDE
    ax1.set_xlabel('Binned correlation')
    ax1.set_ylabel('Density')
    ax1.legend()

    ax2 = fig.add_subplot(122)
    #sns.distplot(open_open_hic, kde=True, hist=True, ax = ax2, label = 'Open - Open (Hi-C)')
    #sns.distplot(open_closed_hic, kde = True, hist=True, ax = ax2, label = 'Open - Closed (Hi-C)')
    #sns.distplot(closed_closed_hic, kde = True, hist=True, ax = ax2, label = 'Closed - Closed (Hi-C)')
    sns.histplot(open_open_hic, kde=True, ax = ax2, label = 'Open - Open (Hi-C)')
    sns.histplot(open_closed_hic, kde = True, ax = ax2, label = 'Open - Closed (Hi-C)')
    sns.histplot(closed_closed_hic, kde = True, ax = ax2, label = 'Closed - Closed (Hi-C)')
    ax2.set_title("Binned correlations (Hi-C)")#Gaussian KDE
    ax2.set_xlabel('Binned correlation')
    ax2.set_ylabel('Density')
    ax2.legend()

    plt.suptitle('Densities of binned correlations', fontsize = 15)
    #plt.show()
    figname = f'kde_binned_corr_chr{args.chrom}_{args.cpg_type}.png'
    fig.tight_layout()
    plt.savefig(os.path.join(args.save_dir, figname))
    print(os.path.join(args.save_dir, figname))
    plt.clf()

