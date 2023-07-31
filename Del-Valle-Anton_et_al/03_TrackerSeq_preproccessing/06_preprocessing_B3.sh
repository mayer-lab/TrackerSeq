#!/bin/bash
#SBATCH -p fat
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -t 48:00:00
#SBATCH --job-name pre3
#SBATCH -o %j.out
#SBATCH -e %j.err


module load bbmap
module load samtools
module load anaconda3
module load openjdk

source ~/.bashrc
conda activate /usr/users/mmnb1ngs/collab_work/workspace/lineage_analysis/conda_lineage


./preprocessing2.sh /home/mpg08/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/fastq_outs/P607_20221213_KBE/3 /usr/users/mmnb1ngs/collab_work/workspace/ykbabal/lin_analysis/corrected_outputs/b3 ED221124_b3_S4_L001_2_R1_001.fastq.gz ED221124_b3_S4_L001_2_R2_001.fastq.gz b3 3000
