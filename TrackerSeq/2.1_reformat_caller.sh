
output=$1
fastq=$2
sample=$3

python 2_reformat_v2.py $output/$fastq $output/barcode_reformat $sample
