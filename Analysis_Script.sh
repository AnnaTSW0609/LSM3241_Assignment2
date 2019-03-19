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
samtools view -f 4 -F 13 results/bam/A0192658N-sorted.bam > results/bam/itself_unmapped_filtered_1.bam 
samtools view -f 8 -F 13 results/bam/A0192658N-sorted.bam > results/bam/mate_unmapped_filtered_2.bam 
# filter both unmapped 
# cannot filter 13 

# command for bedtools to fastq 
# https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html
bedtools bamtofastq -i results/bam/itself_unmapped_filtered.bam -fq results/bam/itself_unmapped_filtered.fq
bedtools bamtofastq -i results/bam/mate_unmapped_filtered.bam -fq results/bam/mate_unmapped_filtered.fq


# check if these can be mapped to the transposome genome 

for file in results/bam/*1.fq ; do
    SRR=$(basename $file _1.fq)
    echo running $SRR
    bowtie2 -x ty5_6p \
         --very-fast -p 4\
         -1 results/bam/itself_unmapped_filtered.fq \  
         -2 results/bam/mate_unmapped_filtered.fq \
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



