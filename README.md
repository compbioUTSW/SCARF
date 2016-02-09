Requirements
------------

Perl - http://www.perl.org
BioPerl - http://www.bioperl.org/wiki/Main_Page
Blat - http://genome.ucsc.edu/FAQ/FAQblat.html#blat3
BWA - http://bio-bwa.sourceforge.net
SAMtools - http://samtools.sourceforge.net

Example commands
----------------

~/scripts/fusion.detection.sh ~/reference/hg19.fasta osteosarcoma {Blood,Primary,Secondary,Tertiary}.realigned.bam
~/scripts/fusion.frequency.sh ~/reference/hg19.fasta osteosarcoma.fusion.txt Secondary 100 Secondary/*.fastq.gz
