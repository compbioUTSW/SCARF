# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use List::Util qw(max sum);
use Bio::DB::Fasta;

my ($positionFile, $referenceFastaFile, @bamFileList) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);
open(my $reader, $positionFile);
while(my $line = <$reader>) {
	chomp($line);
	my ($chromosome, $position, $direction) = split(/\t/, $line);
	my @baseCountHashList = ();
	foreach my $bamFile (@bamFileList) {
		open(my $reader, "samtools view -F 3076 -q 1 $bamFile $chromosome:$position-$position |");
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'} = split(/\t/, $line);
			if($direction eq '-' && $tokenHash{'pos'} == $position && $tokenHash{'cigar'} =~ /^([0-9]+)S/) {
				my $sequence = substr($tokenHash{'seq'}, 0, $1);
				($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
				my @baseList = split(//, $sequence);
				$baseCountHashList[$_]->{$baseList[$_]}++ foreach(0 .. $#baseList);
			}
			if($direction eq '+' && getEnd(@tokenHash{'pos', 'cigar'}) == $position && $tokenHash{'cigar'} =~ /([0-9]+)S$/) {
				my $sequence = substr($tokenHash{'seq'}, -$1);
				my @baseList = split(//, $sequence);
				$baseCountHashList[$_]->{$baseList[$_]}++ foreach(0 .. $#baseList);
			}
		}
		close($reader);
	}
	my @baseList = ();
	my @countList = ();
	my $mismatch = 0;
	foreach my $baseCountHash (@baseCountHashList) {
		my $maximumCount = max(values %$baseCountHash);
		my @maximumCountBaseList = grep {$baseCountHash->{$_} == $maximumCount} keys %$baseCountHash;
		if(scalar(@maximumCountBaseList) == 1) {
			push(@baseList, @maximumCountBaseList);
		} else {
			push(@baseList, 'N');
		}
		push(@countList, sum(map {$baseCountHash->{$_}} @maximumCountBaseList));
		$mismatch += sum(0, map {$baseCountHash->{$_}} grep {$_ ne $maximumCountBaseList[0]} keys %$baseCountHash);
	}
	my @indexList = grep {$countList[$_] > 1} 0 .. $#countList;
	if((my $denominator = sum(0, map {values %$_} @baseCountHashList[@indexList])) > 0) {
		$mismatch = $mismatch / $denominator;
	} else {
		$mismatch = 1;
	}
	my $sequence = '';
	if($direction eq '-') {
		$sequence = uc($db->seq($chromosome, $position, $position + 100 - 1));
		($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
	}
	if($direction eq '+') {
		$sequence = uc($db->seq($chromosome, $position - 100 + 1, $position));
	}
	print join("\t", $chromosome, $position, $direction, scalar(@indexList), $mismatch, $sequence, join('', @baseList), join(',', @countList)), "\n" if(length($sequence) == 100);
}
close($reader);

sub getEnd {
	my ($position, $cigar) = @_;
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		$position += $length if($operation eq 'M');
		$position += $length if($operation eq 'D');
		$position += $length if($operation eq 'N');
	}
	return $position - 1;
}
