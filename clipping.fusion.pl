# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;

my ($pslFile) = @ARGV;
my %filterHash = ();
{
	chomp(my @nameList = `awk -F'\\t' '(NR > 5 && \$12 <= 90 && \$13 >= 110)' $pslFile | cut -f10 | uniq`);
	$filterHash{$_} = 1 foreach(@nameList);
}
my %selectHash = ();
{
	chomp(my @nameList = `awk -F'\\t' '(NR > 5 && \$12 <= 100 && \$13 >= 110)' $pslFile | cut -f10 | uniq -u`);
	$selectHash{$_} = 1 foreach(@nameList);
}
open(my $reader, "awk -F'\\t' '(NR > 5 && \$12 <= 100 && \$13 >= 110 && \$18 == 1)' $pslFile |");
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line);
	unless($filterHash{$tokenList[9]}) {
		if($selectHash{$tokenList[9]}) {
			my $position = 0;
			if($tokenList[8] eq '+') {
				$position = $tokenList[15] + 100 - $tokenList[11] + 1;
			}
			if($tokenList[8] eq '-') {
				$position = $tokenList[16] + $tokenList[11] - 100;
			}
			print join("\t", split(/:/, $tokenList[9]), $tokenList[13], $position, $tokenList[8]), "\n";
		}
	}
}
close($reader);
