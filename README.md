# Inference of tissue type-specific consensus 3D genome using DNA methylation data

## 1. Overview
- Implementation of [Fortin, Jean-Philippe, and Kasper D. Hansen. "Reconstructing A/B compartments as revealed by Hi-C using long-range correlations in epigenetic data." Genome biology 16.1 (2015): 1-23.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0741-y)
- Main idea:
  - Inferring information about 3D genome (A/B compartment distribution) using epigenetic data (450K methylation data)
## 2. Data
- Hi-C data: [GSE18199](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199) 
  - HiC-EBV-2009
  - GSE18199_eigenvector_files.zip.gz
- 450K DNA methylation data: [GSE36369](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36369)
  - GSE36369_series_matrix.txt.gz
