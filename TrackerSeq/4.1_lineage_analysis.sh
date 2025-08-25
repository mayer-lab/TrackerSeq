
dir=$1
sample=$2
n_reads=$3
n_umi=$4
    
python 4_LARRY_Klein.py $dir/${sample}_whitelist.txt $dir/barcode_reformat.fa $dir/${sample}_h5_bartender_pcr_barcode.csv $dir/${sample}_h5_bartender_pcr_cluster.csv  $dir/${sample}_matrix_UMI5 $dir/${sample} $n_reads $n_umi

python 5_cellbc_assign.py $dir/${sample}_matrix_UMI5 $dir/${sample}_whitelist.txt $dir/${sample}_sparse_matrix_UMI5

