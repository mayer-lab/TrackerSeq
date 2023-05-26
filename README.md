# TrackerSeq

TrackerSeq is a 10X Genomics compatible method that labels cells with a unique lineage barcode. This repository contains an updated version of TrackerSeq analysis scripts including multi-core processing that can take raw FASTQ files and convert them into clonal annotations. In this version, Hamming distance calculation was updated with Levenshtein distance calculation in order to increase performance, and cloneID assignment was updated with a network-based approach instead of Jaccard distance calculation. These updates significantly increased the performance of the TranckerSeq pipeline for larger datasets. The original repository of the TrackerSeq pipeline can be found [here](https://github.com/mayer-lab/Bandler-et-al_lineage).

## Dependencies

* [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)
* [Samtools](https://www.htslib.org/)
* [UMI-tools](https://umi-tools.readthedocs.io/en/latest/index.html)
* Python libraries
  * NumPy
  * pandas
  * matplotlib
  * itertools
  * Levenshtein

## TrackerSeq workflow

 The process is divided mainly into four steps: 

1. Pre-processing and barcode matching (**1_pre_process.sh**)
   * Input: Raw Read1 and Read2 files generated from next generation sequencing 
   * Output: 
     1. A fastq file with lineage barcode sequences and headers indicating cell barcode and UMI. 
     2. Cell barcode whitelist .txt file. 
2. Reformat (**2_reformat.ipynb**)
   * Input: The output of step 1 (fastq files of barcodes with headers indicating cell barcode and UMI)
     * ![2_before_reformat](images/2_before_format.png)
   * Output:
     1. Reformmated FASTQ file that's LARRY compatible. 
        * ![2_after_reformat](images/2_after_reformat.png)
     2. A library .txt where each row indicates dataset origin of the lineage barcode. If done correctly, each row should have the same dataset. 
3. Lineage Barcode identification (**3_LARRY_Klein_v2.py**)
   * Input: 
     1. The reformatted FASTQ file that's been converted to a LARRY compatible format from step 2.
     2. Library .txt file generated from step 2.
     3. A .txt file where each row is a cellbc. You can use the whitelist generated in step 1 or get the cell barcodes from your Seurat object.
   * Output: 
     1. CSV file in the form of an NxM binary matrix, where entry (i,j) is 1 if cell i is in clone j
     2. Numpy version of the NxM binary matrix
     3. A barcode .txt file of all the collapsed lineage barcodes
4. Matrix trimming and cellbc assignment (**4_cellbc_assign.ipynb**)
   * Input:
     1. NxM binary matrix CSV file from step 3
     2. Numpy version of the NxM binary from step 3
     3. Cell barcode whitelist generated from step 1
   * Output:
     1. Trimmed matrix CSV file with cellbc column added
5. CloneID assignment (**5_cloneID_assignment.R**)
   * Here we model clonal relationships as a network where nodes are cell-IDs (i.e. cells), which are connected if they share at least one common lineage barcode. The algorithm iterates through all lineage barcodes contained in the sparse matrix from step 4 and progressively updates the edges of the network. Clonal identities are inferred by calculating connected components of the network. This procedure yields equivalent results the the former clustering approach using Jaccard distance (when using maximum distance smaller 1 as a cutoff, ~0.999), but is significantly faster.
   * Input: 
     1. Trimmed sparse matrix CSV file from step 4.
   * Output:
     1. CSV table file with cellbc and corresponding clone ids.


