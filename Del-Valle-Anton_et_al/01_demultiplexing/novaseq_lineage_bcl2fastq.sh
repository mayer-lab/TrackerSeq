#!/bin/bash
#SBATCH -p fat
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -t 18:00:00
#SBATCH --job-name bcl2fastq
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ykbabal@gtu.edu.tr

module load gcc
module load bcl2fastq

bcl2fastq --use-bases-mask Y*,I8n*,I8n*,Y* \
        --minimum-trimmed-read-length=8 \
        --mask-short-adapter-reads=8 \
        --ignore-missing-positions \
        --ignore-missing-controls \
        --ignore-missing-filter \
        --ignore-missing-bcls \
        --no-bgzf-compression \
        --barcode-mismatches 1 \
        --tiles s_1,s_2 \
        -r 15 -w 15 \
        -R /usr/users/mmnb1ngs/collab_work/datasets/sec/P607_221214_A01878_0032_AHLM7GDRX2 \
        --output-dir=/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/fastq_outs/transcriptome \
        --interop-dir=/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/fastq_outs/transcriptome \
        --sample-sheet=/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/SampleSheet_lineage.csv \


