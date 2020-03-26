#!/usr/bin/perl
#For GFF3 files entry

my $fonteDestino = $ARGV[0];
my $gravaDestino = $ARGV[1];
my $line;
open (IN, $fonteDestino) or die "cannot open file";
open(my $fh, '>', $gravaDestino) or die "cannot open file";
while ($line = <IN>){

	chomp($line);
	if ($line=~/.*gene.*Name=(.*)/){
	   print ($fh "$1\n");
	
	}


}
close ($fh);
close (IN);
exit 0;
