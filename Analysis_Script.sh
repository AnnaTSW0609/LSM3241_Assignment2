""" Analysis script used for generating files for Assignement 2"""
""" Raw data files were downloaded as provided from IVLE and were originally in the data folder"""
""" File names had been changed after they are generated for clarity (name changes specified in comments."""

# data preprocessing

# Quality check using FASTQC
mkdir -p results/fastqc
fastqc data/*.fq* -o results/fastqc

# get the text output 
cd results/fastqc
for filename in *.zip; do 
	unzip $filename
done

# document the FASTQC output 
cd ../.. 
mkdir -p docs
cat results/fastqc/*/summary.txt > docs/fastq_summaries.txt

# Trimming bases of poor quality using Trimmomatic

cd data
dpkg -L trimmomatic
cp /usr/share/trimmomatic/NexteraPE-PE.fa . 
                A0192658N_1.trim.fq A0192658N_1un.trim.fq \
                A0192658N_2.trim.fq A0192658N_2un.trim.fq\
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 

gzip A0192658N_1.fq
gzip A0192658N_2.fq


# make a separate directory for the trimmed pairs and orphaned results
mkdir -p results/trimmed
mkdir -p results/orphaned

for file in data/*_1.fq.gz; do
    SRR=$(basename $file _1.fq.gz)
    echo working on $SRR
    Trimmomatic PE data/${SRR}_1.fq.gz data/${SRR}_2.fq.gz \
                  results/trimmed/${SRR}_1.trim.fq.gz results/orphaned/${SRR}_1.untrim.fq.gz \
                  results/trimmed/${SRR}_2.trim.fq.gz results/orphaned/${SRR}_2.untrim.fq.gz \
                  SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
done

# alignment: see if the pair is a proper pair (if each of them maps to a different chr)

# load the genome 

mkdir -p data/genome
mkdir -p results/sam
mkdir -p results/bam
cd data
mv sacCer3.fa genome
mv ty5_6p.fa genome
bowtie2-build data/genome/sacCer3.fasta data/genome/sacCer3 # need to change fa to fasta
bowtie2-build data/genome/ty5_6p.fasta data/genome/ty5_6p


export BOWTIE2_INDEXES=$(pwd)/data/genome

# mapping to the yeast genome 
for file in results/trimmed/*1.trim.fq.gz ; do
    SRR=$(basename $file _1.trim.fq.gz)
    echo running $SRR
    bowtie2 -x sacCer3 \
         --very-fast -p 4\
         -1 results/trimmed/${SRR}_1.trim.fq.gz \
         -2 results/trimmed/${SRR}_2.trim.fq.gz \
         -S results/sam/${SRR}.sam
done
for file in results/sam/*.sam; 
do
	SRR=$(basename $file .sam)
                 echo $SRR
                 samtools view -S -b results/sam/${SRR}.sam > results/bam/${SRR}-aligned.bam
done

for file in results/bam/*-aligned.bam 
do
	SRR=$(basename $file -aligned.bam)
                 echo $SRR
                 samtools sort results/bam/${SRR}-aligned.bam -o results/bam/${SRR}-sorted.bam
done

# check if have one end mapped in yeast genome
samtools view -S -b -f 4 -F 264 results/bam/A0192658N-sorted.bam > results/bam/trial.bam # its mate mapped to the yeast genome, itself unmapped 
samtools view -f 8 -F 260 results/bam/A0192658N-sorted.bam > results/bam/trial_subseq.bam # itself mapped to the yeast genome, its mate unmapped
 
# Merging the two bam files
samtools merge results/bam/merged_unmapped.bam results/bam/trial.bam results/bam/trial_subseq.bam

# convert the singlets back into fasq

samtools sort -n results/bam/merged_unmapped.bam -o results/bam/merged_unmapped.qsort

bedtools bamtofastq -i results/bam/merged_unmapped.qsort \
                      -fq results/trimmed/yeast_unmapped_01.fq \ # the yeast-unmapped are to be mapped to the transposome 
                      -fq2 results/trimmed/yeast_unmapped_02.fq

# mapping the singlets.fastq to the transposon sequence 

for file in results/trimmed/*unmapped_01.fq ; do
    SRR=$(basename $file _unmapped_01.fq)
    echo running $SRR
    bowtie2 -x ty5_6p \
         --very-fast -p 4\
         -1 results/trimmed/${SRR}_unmapped_01.fq \
         -2 results/trimmed/${SRR}_unmapped_02.fq \
         -S results/sam/${SRR}.sam
done

for file in results/sam/*.sam 
do
	SRR=$(basename $file .sam)
                 echo $SRR
                 samtools view -S -b results/sam/${SRR}.sam > results/bam/${SRR}-aligned.bam
done

for file in results/bam/*-aligned.bam 
do
	SRR=$(basename $file -aligned.bam)
                 echo $SRR
                 samtools sort results/bam/${SRR}-aligned.bam -o results/bam/${SRR}-sorted.bam
done

# view the reads that are BOTH singly mapped to the yeast genome and singly mapped to the transposon sequence 
samtools view -S -b -f 4 -F 264 results/bam/yeast-sorted.bam > results/bam/transposome_01.bam # one end mapped on the transposon sequence, itself unmapped  
samtools view -S -b -f 8 -F 260 results/bam/yeast-sorted.bam > results/bam/transposome_02.bam # itself mapped to the transposon sequence, its mate unmapped 

bedtools bamtobed -i results/bam/transposome_02.bam > results/bed/transposome_02.bed

# Merging the two bam files
samtools merge results/bam/merged_trans_unmapped.bam results/bam/transposome_01.bam results/bam/transposome_02.bam # this one are the ones mapped to the transposons 

# convert the singlets back into fasq
samtools sort -n results/bam/merged_trans_unmapped.bam -o results/bam/merged_trans_unmapped.qsort

bedtools bamtofastq -i results/bam/merged_trans_unmapped.qsort \
                      -fq results/trimmed/trans_unmapped_01.fq \ 
                      -fq2 results/trimmed/trans_unmapped_02.fq

# then map these reads to the yeast genome again (running the sam-generating code again)
for file in results/trimmed/*s_unmapped_01.fq ; do
    SRR=$(basename $file ns_unmapped_01.fq)
    echo running $SRR
    bowtie2 -x sacCer3 \
         --very-fast -p 4\
         -1 results/trimmed/${SRR}ns_unmapped_01.fq \
         -2 results/trimmed/${SRR}ns_unmapped_02.fq \
         -S results/sam/${SRR}.sam
done

for file in results/sam/*.sam 
do
	SRR=$(basename $file .sam)
                 echo $SRR
                 samtools view -S -b results/sam/${SRR}.sam > results/bam/${SRR}-aligned.bam
done

for file in results/bam/*-aligned.bam 
do
	SRR=$(basename $file -aligned.bam)
                 echo $SRR
                 samtools sort results/bam/${SRR}-aligned.bam -o results/bam/${SRR}-sorted.bam
done


# obtain for those that are singly mapped to the yeast AND with its mate singly mapped to the transposon
samtools view -S -b -F 4 -F 264 results/bam/tra-sorted.bam > results/bam/tra-sorted_filtered.bam 

# convert it into bed file for creating custom track on the UCSC genome browser  
mkdir results/bed
bedtools bamtobed -i results/bam/tra-sorted_filtered.bam  > results/bed/tra-sorted_filtered.bed #succeed!!! need the headers!!!!!
