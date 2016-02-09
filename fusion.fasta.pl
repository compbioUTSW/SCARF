# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use Bio::DB::Fasta;

my ($fusionFile, $referenceFastaFile, $length) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);
my %nameHash = ();
open(my $reader, $fusionFile);
while(my $line = <$reader>) {
	chomp($line);
	my ($chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2) = split(/\t/, $line);
	my ($name1, $name2, $fusionName) = getNameList($chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2);
	my ($sequence1, $sequence2, $fusionSequence) = getSequenceList($chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2);
	if($fusionSequence ne '') {
		unless(defined($nameHash{$name1})) {
			$nameHash{$name1} = 1;
			print ">$name1\n";
			print "$sequence1\n";
		}
		unless(defined($nameHash{$name2})) {
			$nameHash{$name2} = 1;
			print ">$name2\n";
			print "$sequence2\n";
		}
		print ">$fusionName\n";
		print "$fusionSequence\n";
	}
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

sub getSequenceList {
	my ($chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2) = @_;
	my $sequence1 = '';
	$sequence1 = uc($db->seq($chromosome1, $position1 - ($length - 1), $position1 + $length)) if($strand1 eq '+');
	($sequence1 = reverse(uc($db->seq($chromosome1, $position1 - $length, $position1 + ($length - 1))))) =~ tr/ACGT/TGCA/ if($strand1 eq '-');
	my $sequence2 = '';
	$sequence2 = uc($db->seq($chromosome2, $position2 - $length, $position2 + ($length - 1))) if($strand2 eq '+');
	($sequence2 = reverse(uc($db->seq($chromosome2, $position2 - ($length - 1), $position2 + $length)))) =~ tr/ACGT/TGCA/ if($strand2 eq '-');
	my $fusionSequence = '';
	$fusionSequence = substr($sequence1, 0, $length) . substr($sequence2, $length) if(length($sequence1) == $length * 2 && length($sequence2) == $length * 2);
	return ($sequence1, $sequence2, $fusionSequence);
}
