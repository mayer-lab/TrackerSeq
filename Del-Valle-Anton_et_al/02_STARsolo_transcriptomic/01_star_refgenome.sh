#!/bin/bash
#SBATCH -p fat
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 18:00:00
#SBATCH --job-name star_ref
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ykbabal@gtu.edu.tr




module load star

STAR --runThreadN 20 \
        --runMode genomeGenerate \
        --genomeDir /usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/musFur1_genome_ncbi \
        --genomeFastaFiles /usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/genome/ncbi_genome/ncbi_dataset/data/GCF_000215625.1/GCF_000215625.1_MusPutFur1.0_genomic.fna \
        --sjdbGTFfile /usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/genome/ncbi_genome/Mustelaputorius_enriched.gtf
