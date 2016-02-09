# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/sh

SCRIPT_DIR="`dirname \"$0\"`"
GENOME_FASTA=$1
PREFIX=$2
BAM_FILES=${@:3}

if [ -z "$BAM_FILES" ]; then
	echo 'Usage: fusion.detection.sh <genome.fasta> <output.prefix> <input.bam> ...' 1>&2
	exit 1
fi
perl -MBio::DB::Fasta -e '' || exit 1
type blat > /dev/null || exit 1

perl $SCRIPT_DIR/clipping.position.pl $BAM_FILES > $PREFIX.clipping.position.txt

# Clipping reads >= 2
awk -F'\t' '($4 >= 2)' $PREFIX.clipping.position.txt | perl $SCRIPT_DIR/clipping.sequence.pl - $GENOME_FASTA $BAM_FILES > $PREFIX.clipping.sequence.txt

# Length covered by at least 2 reads >= 10 and mismatch ratio in the length < 0.1
awk -F'\t' '($4 >= 10 && $5 < 0.1)' $PREFIX.clipping.sequence.txt | awk -F'\t' -vOFS='\n' '{print ">"$1":"$2":"$3, $6, $7}' > $PREFIX.clipping.fasta

blat $GENOME_FASTA $PREFIX.clipping.fasta $PREFIX.clipping.blat.psl

perl $SCRIPT_DIR/clipping.fusion.pl $PREFIX.clipping.blat.psl | perl $SCRIPT_DIR/fusion.unique.pl - $GENOME_FASTA > $PREFIX.fusion.txt
