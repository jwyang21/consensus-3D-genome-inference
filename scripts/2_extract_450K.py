import pandas as pd
import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', required=True, default = '14', help = 'chromosome-of-interest', type = str)
    parser.add_argument('-w_dir', '--working_dir', required = True, help = 'working directory', default = '/data/project/jeewon/research/3dith-reproduce', type = str)
    parser.add_argument('-d_dir', '--data_dir', required = True, help = '450K data directory', default = '/data/project/jeewon/research/3dith-reproduce/data', type = str)
    parser.add_argument('-s_dir', '--save_dir', required = True, help = 'saving directory', default = '/data/project/jeewon/research/3dith-reproduce/result', type = str)
    parser.add_argument('--cpg_type', required = True, help = 'island, opensea, shelf, shore', type = str)
    parser.add_argument('--grch', required = True, help = 'GRCh version', default = 36, type = int)
    return parser.parse_args()

suppl_fname = 'GPL13534_HumanMethylation450_15017482_v.1.1.csv'
manifest_fname = 'GPL13534_450K_Manifest_header_Descriptions.xlsx'
data_fname = 'GSE36369_series_matrix.txt'

if __name__ == '__main__':
    
    args = parse_arguments()
    if not os.path.exists(args.data_dir):
        os.makedirs(args.data_dir)
    if not os.path.exists(args.working_dir):
        os.makedirs(args.working_dir) 
    
    os.chdir(args.working_dir)
    
    if data_fname not in os.listdir(args.data_dir):
        os.chdir(args.data_dir)
        cmd1_ = 'wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE36nnn/GSE36369/matrix/GSE36369_series_matrix.txt.gz'
        cmd2_ = 'gzip -d GSE36369_series_matrix.txt.gz'
        os.system(cmd1_)
        os.system(cmd2_)
        os.chdir(args.working_dir)
    
    # 1. Import 450K beta value
    ## 1-1. Load beta value
    data_fname = os.path.join(args.data_dir, data_fname)
    f = open(data_fname, 'r')
    f_content = []
    while True:
        line = f.readline()
        if not line:
            break
        if line.startswith('!')==False:
            if line.startswith("\n")==False:
                tmp = line.split('\t')
                tmp2 = [x.strip() for x in tmp]
                f_content.append(tmp2)
    f.close()
    f_content_arr = np.array(f_content)
    data = pd.DataFrame(f_content_arr[1:,1:], index = f_content_arr[1:,0].flatten(), columns = f_content_arr[0,1:].flatten())
    print('---\n450K beta value data: ')
    print(data.head(3))
    print(data.shape)
    
    ## 1-2. convert data type (str -> float)
    data.replace('', np.nan, inplace=True)
    data = data.astype('float')
    
    ## 1-3. arrange index of data
    new_index = [x.strip('"') for x in data.index.values]
    data.index = new_index
    
    ## 1-4. save beta value
    result_fname = f'450K-GRCh{args.grch}-beta.csv'
    data.to_csv(os.path.join(args.save_dir, result_fname))
    print(os.path.join(args.save_dir, result_fname))
    
    # 2. Import target cpg probe list
    cpg_fname = os.path.join(args.save_dir, f'450K-GRCh{args.grch}-chr{args.chrom}-{args.cpg_type}.csv')
    cpg = pd.read_csv(cpg_fname, index_col = 0)
    
    # 3. extract beta value of target cpg probes
    data2 = data[data.index.isin(cpg.index.values)].copy()
    data2_fname = os.path.join(args.save_dir, f'450K-GRCh{args.grch}-chr{args.chrom}-{args.cpg_type}-beta.csv')
    data2.to_csv(data2_fname)
    print(data2_fname)
