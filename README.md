Requirements
------------

Perl - http://www.perl.org

BioPerl - http://www.bioperl.org/wiki/Main_Page

Blat - http://genome.ucsc.edu/FAQ/FAQblat.html#blat3

BWA - http://bio-bwa.sourceforge.net

SAMtools - http://samtools.sourceforge.net

Example commands
----------------

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
