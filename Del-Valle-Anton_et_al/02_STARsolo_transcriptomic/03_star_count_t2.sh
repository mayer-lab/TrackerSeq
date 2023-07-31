#!/bin/bash
#SBATCH -p fat
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 18:00:00
#SBATCH --job-name star_t2
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ykbabal@gtu.edu.tr




module load anaconda3
source ~/.bashrc
conda activate /usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/conda_star


fastq=/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/fastq_outs/transcriptome/P607_20221213_KBE/2

STAR --runThreadN 20 \
        --genomeDir /usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/musFur1_genome_ncbi \
        --readFilesIn $fastq/ED221124_t2_S3_L001_R2_001.fastq.gz,$fastq/ED221124_t2_S3_L002_R2_001.fastq.gz $fastq/ED221124_t2_S3_L001_R1_001.fastq.gz,$fastq/ED221124_t2_S3_L002_R1_001.fastq.gz \
        --soloType CB_UMI_Simple \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 12 \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes CR UR CY UY CB UB \
        --soloCBwhitelist /usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/res_STAR/3M-february-2018.txt \
        --readFilesCommand gunzip -c
