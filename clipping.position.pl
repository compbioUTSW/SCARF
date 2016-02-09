# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use Getopt::Long;

GetOptions('r=s' => \(my $region = ''));
my (@bamFileList) = @ARGV;
open(my $writer, "| sort --field-separator=\$'\\t' -k1,1 -k2,2n -k3,3 | uniq -c | sed 's/^ *\\([0-9]*\\) \\(.*\\)\$/\\2\\t\\1/'");
foreach my $bamFile (@bamFileList) {
	open(my $reader, "samtools view -F 3076 -q 1 $bamFile $region |");
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'} = split(/\t/, $line);
		print $writer join("\t", $tokenHash{'rname'},        $tokenHash{'pos'},           '-'), "\n" if($tokenHash{'cigar'} =~ /^[0-9]+S/);
		print $writer join("\t", $tokenHash{'rname'}, getEnd(@tokenHash{'pos', 'cigar'}), '+'), "\n" if($tokenHash{'cigar'} =~ /[0-9]+S$/);
	}
	close($reader);
}
close($writer);

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
