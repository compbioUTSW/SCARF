Description
-----------
This folder contains scripts that can be used to detect structural variants such as inter- or intra-chromosomal translocations, which may generate abnormal fusion transcripts, from exome-sequencing data and esitmate their frequency.  

A pictorial description of the workflow can be found in [workflow.pdf](https://github.com/compbioUTSW/fusion/blob/master/workflow.pdf), and a synthetically generated toy example is also provided as [BCR_ABL1.fastq.gz](https://github.com/compbioUTSW/fusion/blob/master/example/BCR_ABL1.fastq.gz).  

For further questions please contact [@jiwoongbio](https://github.com/jiwoongbio).

Requirements
------------

Perl - http://www.perl.org

BioPerl - http://www.bioperl.org/wiki/Main_Page

Blat - http://genome.ucsc.edu/FAQ/FAQblat.html#blat3

BWA - http://bio-bwa.sourceforge.net

SAMtools - http://samtools.sourceforge.net

Example commands
----------------
```
# Move to directory "example"
cd example

# Download human genome sequence and build the index
# - You can replace this step with using your prebuilt human genome index.
wget ftp-trace.ncbi.nih.gov:/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
gzip -d human_g1k_v37.fasta.gz
bwa index human_g1k_v37.fasta

# Mapping sequencing reads to human genome
bwa mem human_g1k_v37.fasta BCR_ABL1.fastq.gz | samtools view -S -b - > BCR_ABL1.bam

# Sort alignment file and build the index
samtools sort BCR_ABL1.bam BCR_ABL1.sorted
samtools index BCR_ABL1.sorted.bam

# Fusion detection
../fusion.detection.sh human_g1k_v37.fasta BCR_ABL1 BCR_ABL1.sorted.bam

# Fusion frequency estimation
../fusion.frequency.sh human_g1k_v37.fasta BCR_ABL1.fusion.txt BCR_ABL1 100 BCR_ABL1.fastq.gz
```

Citation
----------------
A novel TP53-KPNA3 translocation defines a de novo treatment-resistant clone in osteosarcoma.
Cold Spring Harbor molecular case studies. 2016 Sep;2(5):a000992.
http://molecularcasestudies.cshlp.org/content/2/5/a000992.full
