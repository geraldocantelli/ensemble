#!/usr/bin/perl
#For GTF files entry

my $fonteDestino = $ARGV[0];
my $gravaDestino = $ARGV[1];
my $line;
open (IN, $fonteDestino) or die "cannot open file";
open(my $fh, '>', $gravaDestino) or die "cannot open file";
while ($line = <IN>){

	chomp($line);
	if ($line=~/.*gene_name "([\w\d\_]+)";(.*)/){
	   print ($fh "$1\n");
	
	}


}
close ($fh);
close (IN);
exit 0;

