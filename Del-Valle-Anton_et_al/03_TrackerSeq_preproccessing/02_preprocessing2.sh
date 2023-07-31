datadir=$1
output=$2
R1=$3
R2=$4
sample=$5
expect=$6


### 1.Trimming and selecting
# Replace in1 and in2 with your R1 and R2 of your retrieved barcode file
bbduk.sh in=$datadir/$R1 in2=$datadir/$R2 k=17 literal=GACTCTGGCTCACAAAT  ktrim=r out=stdout.fq int=f skipr1 maq=15 |\
        bbduk.sh in=stdin.fq literal=CTGA k=4 restrictleft=4 ktrim=l out1=$output/trim_R1.fastq out2=$output/trim_R2.fastq outm1=$output/dis_R1.fastq outm2=$output/dis_R2.fastq int=t skipr1

### 2. Remove short reads and re-pair
# Remove short reads and re-pair R1 and R2 (for libraries using piggybac)

bbduk.sh in=$output/trim_R2.fastq out=stdout.fq minlength=37 maxlength=37 | \
        repair.sh in1=stdin.fq in2=$output/trim_R1.fastq out1=$output/pair_R1.fastq out2=$output/pair_R2.fastq repair


### 3. identify correct cell barcodes

# You may have to adjust the settings for umi_tools whitelist depending on the file for more info see: 
umi_tools whitelist --stdin $datadir/$R1 \
        --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
        --knee-method=density \
        --set-cell-number=$expect \
        --plot-prefix=${sample}_expect_whitelist \
        --log2stderr > $output/${sample}_whitelist.txt

### 4. Add the cell barcodes from step 3 to R2 reads 

umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
        --stdin $output/pair_R1.fastq \
        --stdout $output/pair_R1_extracted.fastq.gz \
        --read2-in $output/pair_R2.fastq \
        --read2-out=$output/barcode_extracted.fastq \
        --filter-cell-barcode \
        --error-correct-cell \
        --whitelist=$output/${sample}_whitelist.txt

