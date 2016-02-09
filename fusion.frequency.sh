# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/sh

SCRIPT_DIR="`dirname \"$0\"`"
GENOME_FASTA=$1
FUSION_FILE=$2
PREFIX=$3
READ_LENGTH=$4
FASTQ_FILES=${@:5}

if [ -z "$FASTQ_FILES" ]; then
	echo 'Usage: fusion.frequency.sh <genome.fasta> <input.fusion> <output.prefix> <read.length> <input.fastq> ...' 1>&2
	exit 1
fi
perl -MBio::DB::Fasta -e '' || exit 1
type bwa > /dev/null || exit 1
type samtools > /dev/null || exit 1

perl $SCRIPT_DIR/fusion.fasta.pl $FUSION_FILE $GENOME_FASTA $READ_LENGTH > $PREFIX.fusion.fasta

bwa index -a is $PREFIX.fusion.fasta

for file in $FASTQ_FILES; do bwa aln $PREFIX.fusion.fasta $file | bwa samse -n 30 $PREFIX.fusion.fasta - $file | samtools view -S -F 4 - | perl $SCRIPT_DIR/bwa/xa2multi.pl; done | gzip > $PREFIX.fusion.sam.gz

perl $SCRIPT_DIR/fusion.count.pl $FUSION_FILE $PREFIX.fusion.sam.gz $READ_LENGTH | awk -F'\t' -vOFS='\t' '($7 + $8 + $9 > 0) {print $1, $2, $3, $4, $5, $6, $7 + $8 + $9, $9 / ($7 + $8 + $9)}' > $PREFIX.fusion.frequency.txt
