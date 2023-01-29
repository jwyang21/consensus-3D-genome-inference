chrom=14
w_dir=/data/project/jeewon/research/3dith-reproduce
d_dir=/data/project/jeewon/research/3dith-reproduce/data
s_dir=/data/project/jeewon/research/3dith-reproduce/result
cpg_type=opensea
grch=36
chrom_len=/data/project/jeewon/research/reference/hg18.fa.sizes
binsize=100000

echo "python3 4_binning_correlation.py --chrom $chrom -w_dir $w_dir -d_dir $d_dir -s_dir $s_dir --cpg_type $cpg_type --grch $grch --chrom_len $chrom_len --binsize $binsize"
python3 4_binning_correlation.py --chrom $chrom -w_dir $w_dir -d_dir $d_dir -s_dir $s_dir --cpg_type $cpg_type --grch $grch --chrom_len $chrom_len --binsize $binsize
