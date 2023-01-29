chrom=14
w_dir=/data/project/jeewon/research/3dith-reproduce
d_dir=/data/project/jeewon/research/3dith-reproduce/data
s_dir=/data/project/jeewon/research/3dith-reproduce/result
cpg_type=opensea
grch=36

echo "python3 1_list_chr-of-interest_opensea_450K.py --chrom $chrom -w_dir $w_dir -d_dir $d_dir -s_dir $s_dir --cpg_type $cpg_type --grch $grch"
python3 1_list_chr-of-interest_opensea_450K.py --chrom $chrom -w_dir $w_dir -d_dir $d_dir -s_dir $s_dir --cpg_type $cpg_type --grch $grch
