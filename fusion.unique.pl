# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use Bio::DB::Fasta;

my ($fusionFile, $referenceFastaFile) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);
chomp(my @chromosomeList = `grep '^>' $referenceFastaFile | sed 's/^>//'`);
s/ .*$// foreach(@chromosomeList);
my %chromosomeIndexHash = map {$chromosomeList[$_] => $_} 0 .. $#chromosomeList;
open(my $reader, $fusionFile);
open(my $writer, "| sort --field-separator=\$'\\t' -k1,1n -k3,3n -k4,7 | cut -f2,3,4,5,6,7 | uniq");
while(my $line = <$reader>) {
	chomp($line);
	my ($chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2) = split(/\t/, $line);
	if(($chromosome1 eq $chromosome2 && $position1 > $position2) || $chromosomeIndexHash{$chromosome1} > $chromosomeIndexHash{$chromosome2}) {
		($chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2) = ($chromosome2, $position2, $strand2, $chromosome1, $position1, $strand1);
		if($strand1 eq '+') {
			$strand1 = '-';
		} elsif($strand1 = '-') {
			$strand1 = '+';
		}
		if($strand2 eq '+') {
			$strand2 = '-';
		} elsif($strand2 = '-') {
			$strand2 = '+';
		}
	}
	($position1, $position2) = leftalignFusion($chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2);
	print $writer join("\t", $chromosomeIndexHash{$chromosome1}, $chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2), "\n";
}
close($reader);
close($writer);

sub leftalignFusion {
	my ($chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2) = @_;
	while(1) {
		my $newPosition1 = $position1 - 1;
		my $newPosition2 = $strand1 eq $strand2 ? $position2 - 1 : $position2 + 1;
		if($strand1 eq '+') {
			my $base1 = uc($db->seq($chromosome1, $position1, $position1));
			my $base2 = uc($db->seq($chromosome2, $newPosition2, $newPosition2));
			$base2 =~ tr/ACGT/TGCA/ if($strand2 eq '-');
			last if($base1 ne $base2);
		} else {
			my $base1 = uc($db->seq($chromosome1, $newPosition1, $newPosition1));
			my $base2 = uc($db->seq($chromosome2, $position2, $position2));
			$base1 =~ tr/ACGT/TGCA/;
			$base2 =~ tr/ACGT/TGCA/ if($strand2 eq '-');
			last if($base1 ne $base2);
		}
		($position1, $position2) = ($newPosition1, $newPosition2);
	}
	return ($position1, $position2);
}
