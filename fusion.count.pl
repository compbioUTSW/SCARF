# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;

my ($fusionFile, $samFile, $length) = @ARGV;
my %nameCountHash = getNameCountHash();
open(my $reader, $fusionFile);
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, ($line =~ tr/\t//) + 1);
	my ($chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2) = @tokenList;
	my ($name1, $name2, $fusionName) = getNameList($chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2);
	print join("\t", @tokenList, map {defined($_) ? $_ : 0} @nameCountHash{$name1, $name2, $fusionName}), "\n";
}
close($reader);

sub getNameList {
	my ($chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2) = @_;
	my $name1 = '';
	$name1 = join(':', $chromosome1, $position1, $strand1, $chromosome1, $position1 + 1, $strand1) if($strand1 eq '+');
	$name1 = join(':', $chromosome1, $position1, $strand1, $chromosome1, $position1 - 1, $strand1) if($strand1 eq '-');
	my $name2 = '';
	$name2 = join(':', $chromosome2, $position2 - 1, $strand2, $chromosome2, $position2, $strand2) if($strand2 eq '+');
	$name2 = join(':', $chromosome2, $position2 + 1, $strand2, $chromosome2, $position2, $strand2) if($strand2 eq '-');
	my $fusionName = join(':', $chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2);
	return ($name1, $name2, $fusionName);
}

sub getNameCountHash {
	my %nameCountHash = ();
	open(my $reader, ($samFile =~ /\.gz$/) ? "gzip -dc $samFile |" : $samFile);
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^@/);
		my %tokenHash = ();
		@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'} = split(/\t/, $line);
		$nameCountHash{$tokenHash{'rname'}}++ if($tokenHash{'pos'} <= $length - 10 + 1 && getEnd(@tokenHash{'pos', 'cigar'}) >= $length + 10);
	}
	close($reader);
	return %nameCountHash;
}

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
