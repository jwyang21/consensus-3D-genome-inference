# Inference of tissue type-specific consensus 3D genome using DNA methylation data

## 1. Overview
- Python implementation of [Fortin, Jean-Philippe, and Kasper D. Hansen. "Reconstructing A/B compartments as revealed by Hi-C using long-range correlations in epigenetic data." Genome biology 16.1 (2015): 1-23.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0741-y)
- Main idea:
  - Inferring information about 3D genome (A/B compartment distribution) using epigenetic data (450K methylation data)
## 2. Data
- 450K DNA methylation data: [GSE36369](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36369)
  - GSE36369_series_matrix.txt.gz
  - Used to infer 3D genome
- Hi-C data: [GSE18199](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199) 
  - HiC-EBV-2009
  - GSE18199_binned_heatmaps.zip.gz
  - Not used in 3D genome inference.
  - Used for comparison only.
## 3. Installation
```shell
conda env create -f consensus-3d-450k.yaml
```
## 4. Run
- Making log files is not mandatory.  
- Commands below are for chromosome 14 and open sea CpG probes. 
  - Values of variables 'chrom' and 'cpg\_type' in shell scripts can be edited to analyze other chromosome and/or other CpG probes.
```shell
cd scripts
bash 1_list_chr-of-interest_opensea_450K.sh > ../log/1_list_chr-of-interest_opensea_450K.log           
bash 2_extract_450K.sh > ../log/2_extract_450K.log         
bash 3_compute_correlation.sh > ../log/3_compute_correlation.log       
bash 4_binning_correlation.sh > ../log/4_binning_correlation.log        
bash 5_compare_hic_450k.sh > ../log/5_compare_hic_450k.log       
```
