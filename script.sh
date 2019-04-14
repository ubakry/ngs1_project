###########################################
#  - NGS Course - Assignment              #
#  - Bash Script                          #
#  - April 1,2019                         #
#  - Copyright: Usama Bakry               #
#  - Nile University                      #
###########################################
#!/bin/bash

# Getting the arguments:
# Input: SRA file
# Output: output directory
while getopts i:o: option
do
case "${option}"
in
i) INPUT=${OPTARG};;
o) OUTPUT=${OPTARG};;
esac
done

# Downloading database from gencode website
# -p => no error if existing, make parent directories as needed
mkdir -p $OUTPUT/db/
# -P => save files to PREFIX/.. (specifying directory for download)
wget -P $OUTPUT/db/ ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh38.p12.genome.fa.gz
wget -P $OUTPUT/db/ ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz


# Printing input and output
echo [INFO] Input:$INPUT
echo [INFO] Output:$OUTPUT

# User confirmation for continue
# y => continue
# n => bash script will stop
# * => invalid
read -p "Continue (y/n)?" choice
case "$choice" in 
  y|Y ) 

# Prepare the data ######################################################################

# Making output directory
echo "[PROCESS] Creating results directory..."
mkdir -p $OUTPUT/results/fastq/

# Converting sra file to R1 and R2 fastq files
echo "[PROCESS] Converting SRA file to R1 and R2 Fastq files..."
# --split-3 => splitting sra file to files *_1.fastq and *_2.fastq for R1 and R2
fastq-dump --split-3 $INPUT -O $OUTPUT/results/fastq
echo "[INFO] Files are converted in "$OUTPUT"/results/fastq/"

# Creating samples (R1 and R2) before shuffling
echo "[PROCESS] Creating samples before shuffling..."
# head => Getting the first 5M reads
# split2 => Splitting the 5M reads into 5 files with 1M read for each one
cat $OUTPUT/results/fastq/$(echo $INPUT | egrep -o SRR[0-9]+)_1.fastq | seqkit head -n 5000000 | seqkit split2 -s 1000000 -O $OUTPUT/results/before/R1
cat $OUTPUT/results/fastq/$(echo $INPUT | egrep -o SRR[0-9]+)_2.fastq | seqkit head -n 5000000 | seqkit split2 -s 1000000 -O $OUTPUT/results/before/R2
echo "[INFO] Samples before shuffling are created in "$OUTPUT"/results/before/"

# Shuffling R1 and R2 fastq files
echo "[PROCESS] Shuffling R1 and R2 fastq files..."
seqkit shuffle --threads 4 $OUTPUT/results/fastq/$(echo $INPUT | egrep -o SRR[0-9]+)_1.fastq > $OUTPUT/results/fastq/$(echo $INPUT | egrep -o SRR[0-9]+)_1_shuffled.fastq
seqkit shuffle --threads 4 $OUTPUT/results/fastq/$(echo $INPUT | egrep -o SRR[0-9]+)_2.fastq > $OUTPUT/results/fastq/$(echo $INPUT | egrep -o SRR[0-9]+)_2_shuffled.fastq
echo "[INFO] Shuffling R1 and R2 is done"

# Creating samples (R1 and R2) after shuffling
echo "[PROCESS] Creating samples after shuffling..."
# head => Getting the first 5M shuffled reads
# split2 => Splitting the 5M shuffled reads into 5 files with 1M read for each one
cat $OUTPUT/results/fastq/$(echo $INPUT | egrep -o SRR[0-9]+)_1_shuffled.fastq | seqkit head -n 5000000 | seqkit split2 -s 1000000 -O $OUTPUT/results/after/R1
cat $OUTPUT/results/fastq/$(echo $INPUT | egrep -o SRR[0-9]+)_2_shuffled.fastq | seqkit head -n 5000000 | seqkit split2 -s 1000000 -O $OUTPUT/results/after/R2
echo "[INFO] Samples after shuffling are created in "$OUTPUT"/results/after/"

# Moving and renaming samples from before and after
echo "[PROCESS] Moving and renaming samples from before and after..."
mkdir -p $OUTPUT/results/samples/before/
mkdir -p $OUTPUT/results/samples/after/
for dir in before after; do
  for i in {1..2}; do
    for j in {1..5}; do
        mv $OUTPUT/results/$dir/R$i/stdin.part_00$j.fastq $OUTPUT/results/samples/$dir/S$j-R$i.fastq
    done
  done
done
echo "[INFO] Samples are moved in "$OUTPUT"/results/samples/"

# -d => remove empty directories
# -r => remove directories and their contents recursively
rm -d -r $OUTPUT/results/before/
rm -d -r $OUTPUT/results/after/

# Prepare the data ######################################################################

# FASTQ Quality Control #################################################################

# Applying quality check on samples using FastQC
echo "[PROCESS] Applying quality check on samples using FastQC..."
mkdir -p $OUTPUT/results/fastqc-reports/before/
mkdir -p $OUTPUT/results/fastqc-reports/after/
# --quiet => Supress all progress messages on stdout and only report errors.
for dir in before after; do
  fastqc $OUTPUT/results/samples/$dir/*.fastq --quiet -o $OUTPUT/results/fastqc-reports/$dir/ 
done
echo "[INFO] Check the quality reports on "$OUTPUT"/results/fastqc-reports/"

# FASTQ Quality Control #################################################################

# Trimming ##############################################################################

# Applying mild trimming on un-shuffled samples using Trimmomatic
echo "[PROCESS] Applying mild trimming on un-shuffled samples using Trimmomatic..."
for j in {1..5}; do
  mkdir -p $OUTPUT/results/trim-files/before/S$j/
  dir=$OUTPUT/results/trim-files/before/S$j
  adapt=/mnt/sda3/Software/trimmomatic/Trimmomatic-0.38/adapters
  trimmomatic PE -phred33 $OUTPUT/results/samples/before/S$j-R1.fastq $OUTPUT/results/samples/before/S$j-R2.fastq $dir/output_R1_paired.fq.gz $dir/output_R1_unpaired.fq.gz $dir/output_R2_paired.fq.gz $dir/output_R2_unpaired.fq.gz ILLUMINACLIP:$adapt/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
echo "[INFO] Mild trimming is done."

# Applying aggressive trimming on un-shuffled samples using Trimmomatic
echo "[PROCESS] Applying aggressive trimming on shuffled samples using Trimmomatic..."
for j in {1..5}; do
  mkdir -p $OUTPUT/results/trim-files/after/S$j/
  dir=$OUTPUT/results/trim-files/after/S$j
  adapt=/mnt/sda3/Software/trimmomatic/Trimmomatic-0.38/adapters
  trimmomatic PE -phred33 $OUTPUT/results/samples/after/S$j-R1.fastq $OUTPUT/results/samples/after/S$j-R2.fastq $dir/output_R1_paired.fq.gz $dir/output_R1_unpaired.fq.gz $dir/output_R2_paired.fq.gz $dir/output_R2_unpaired.fq.gz ILLUMINACLIP:$adapt/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:30 MINLEN:90
done
echo "[INFO] Aggressive trimming is done."

# Trimming ##############################################################################

# Alignment #############################################################################

# Indexing the reference genome
echo "[PROCESS] Indexing the reference genome..."
gunzip -k $OUTPUT/db/GRCh38.p12.genome.fa.gz
mkdir -p $OUTPUT/db/bwa-index/
cd $OUTPUT/db/bwa-index/
ln -s ../GRCh38.p12.genome.fa .
bwa index -a bwtsw GRCh38.p12.genome.fa
cd -
echo "[INFO] Indexing is done"

# Aligning un-shuffled samples using BWA
echo "[PROCESS] Aligning un-shuffled samples using BWA..."
for j in {1..5}; do
  mkdir -p $OUTPUT/results/align-res/bwa/S$j/
  dir=$OUTPUT/results/align-res/bwa/S$j/
  bwa mem $OUTPUT/db/bwa-index/GRCh38.p12.genome.fa $OUTPUT/results/trim-files/before/S$j/output_R1_paired.fq.gz $OUTPUT/results/trim-files/before/S$j/output_R2_paired.fq.gz > $dir/S$j.sam
  samtools flagstat $dir/S$j.sam > $dir/S$j.stat
done
echo "[INFO] BWA Aligning is done."

# Indexing the reference genome
echo "[PROCESS] Indexing the reference genome..."
mkdir -p $OUTPUT/db/hisat-index/
#gunzip -k $OUTPUT/db/gencode.v30.annotation.gtf.gz
cd $OUTPUT/db/hisat-index/
ln -s ../GRCh38.p12.genome.fa .
ln -s ../gencode.v30.annotation.gtf .
hisat2_extract_splice_sites.py gencode.v30.annotation.gtf > splicesites.tsv
hisat2_extract_exons.py gencode.v30.annotation.gtf > exons.tsv
hisat2-build -p 8 --ss splicesites.tsv --exon exons.tsv GRCh38.p12.genome.fa GRCh38.p12.genome
cd -
echo "[INFO] Indexing is done"

# Aligning shuffled samples using Hisat2
echo "[PROCESS] Aligning shuffled samples using Hisat2..."
for j in {1..5}; do
  mkdir -p $OUTPUT/results/align-res/hisat/S$j/
  dir=$OUTPUT/results/align-res/hisat/S$j/
  hisat2 -p 2 -x $OUTPUT/db/hisat-index/GRCh38.p12.genome --dta --rna-strandness RF -1 $OUTPUT/results/trim-files/after/S$j/output_R1_paired.fq.gz -2 $OUTPUT/results/trim-files/after/S$j/output_R2_paired.fq.gz -S $dir/S$j.sam
  samtools flagstat $dir/S$j.sam > $dir/S$j.stat
done
echo "[INFO] Hisat Aligning is done."

# Alignment #############################################################################

# Assembly ##############################################################################

# Preparing the SAM file for assembly
echo "[PROCESS] Preparing the SAM file for assembly..."
for i in bwa hisat; do
  for j in {1..5}; do
      dir=$OUTPUT/results/align-res/$i/S$j/
      # converting the SAM file into BAM file 
      samtools view -bS $dir/S$j.sam > $dir/S$j.bam
      #converting the BAM file to a sorted BAM file. 
      samtools sort $dir/S$j.bam -o $dir/S$j.sorted.bam
  done
done
echo "[INFO] BAM files are located in "$OUTPUT"/results/align-res/"

# Applying assembly on shuffeled and unshuffeled samples
echo "[PROCESS] Applying assembly on shuffeled and unshuffeled samples using StringTie..."
for i in bwa hisat; do
  for j in {1..5}; do
    mkdir -p $OUTPUT/results/assembly-res/$i/S$j/
    dir1=$OUTPUT/results/assembly-res/$i/S$j/
    dir2=$OUTPUT/results/align-res/$i/S$j/
    stringtie $dir2/S$j.sorted.bam -G $OUTPUT/db/gencode.v30.annotation.gtf -o $dir1/S$j.gtf 
  done
done
echo "[INFO] Assembly is done."

# Assembly ##############################################################################

# GTF Compare ###########################################################################

# Comparing GTF files using GFFCompare
echo "[PROCESS] Comparing GTF files using GFFCompare..."
for i in bwa hisat; do
  for j in {1..5}; do
    mkdir -p $OUTPUT/results/compare-res/$i/S$j/
    dir=$OUTPUT/results/compare-res/$i/S$j/
    cd $dir
    gffcompare -r ../../../../db/gencode.v30.annotation.gtf ../../../assembly-res/$i/S$j/S$j.gtf 
    cd -
  done
done
echo "[INFO] Comparing is done."

# GTF Compare ###########################################################################

# Diff. Expression ######################################################################
echo "[PROCESS] Copying files..."
for i in bwa hisat; do
  for j in {1..5}; do
    mkdir -p $OUTPUT/results/diff-res/$i/
    dir1=$OUTPUT/results/align-res/$i/S$j/
    dir2=$OUTPUT/results/diff-res/$i/
    cp $dir1/S$j.bam $dir2/S$j.bam 
  done
done

# Generate the counts.
echo "[PROCESS] Generating the genes counts..."
featureCounts -a $OUTPUT/db/gencode.v30.annotation.gtf -g gene_name -o $OUTPUT/results/diff-res/counts.txt  $OUTPUT/results/diff-res/bwa/*.bam  $OUTPUT/results/diff-res/hisat/*.bam
Simplify the file to keep only the count columns.
cat $OUTPUT/results/diff-res/counts.txt | cut -f 1,7-12 > $OUTPUT/results/diff-res/simple_counts.txt
cat $OUTPUT/results/diff-res/simple_counts.txt | Rscript deseq1.r 5x5 > $OUTPUT/results/diff-res/results_deseq1.tsv
cat $OUTPUT/results/diff-res/results_deseq1.tsv | awk ' $8 < 0.05 { print $0 }' > $OUTPUT/results/diff-res/filtered_results_deseq1.tsv
cat $OUTPUT/results/diff-res/filtered_results_deseq1.tsv | Rscript draw-heatmap.r > $OUTPUT/results/diff-res/hisat_output.pdf

# # Diff. Expression ######################################################################

;;

  n|N ) echo "Process stoped.";;
  * ) echo "invalid";;
esac
