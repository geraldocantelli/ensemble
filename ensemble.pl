use warnings;
use Carp;
use strict;
 
use Venn::Chart;
use Array::Contains;
use List::Uniq ':all';
use Bio::Seq;
use Bio::SeqIO;


#!/usr/bin/perl
my $gravaDestinoSW1 = $ARGV[1];
my $gravaDestinoSW2 = $ARGV[3];
my $INTERSEC = "/home/geraldo/Documentos/Fatec/MestradoUTFPR/src/mestradoutfpr/INTERSEC.txt";

my $onlySW1Perl = "/home/geraldo/Documentos/Fatec/MestradoUTFPR/src/mestradoutfpr/onlySW1Perl.txt";
my $onlySW2Perl = "/home/geraldo/Documentos/Fatec/MestradoUTFPR/src/mestradoutfpr/onlySW2Perl.txt";

my $line;
open (IN, $gravaDestinoSW1) or die "cannot open file";

open(my $fh, '>', $INTERSEC) or die "cannot open file";

my @SW1;
my @SW2;

##################################################################################
#Creating SW1 array

while ($line = <IN>){

	chomp($line);
	if ($line=~/(.*)/){
	   print ("$1\n");
	   push(@SW1,$1);
	}

}

close (IN);

###############################################################################
#Creating SW2 array

open (IN2, $gravaDestinoSW2) or die "cannot open file";

while ($line = <IN2>){

	chomp($line);
	if ($line=~/(.*)/){
	   print ("$1\n");
	   push(@SW2,$1);
	}

}

close (IN2);

################################################################################
#Creating Union set
my @INTERSECTION;
my @uniqSW1= uniq(@SW1);
my @uniqSW2= uniq(@SW2);

foreach my $p (@uniqSW1){

	if(contains($p,\@uniqSW2)){
		push(@INTERSECTION,$p);
		print($fh "$p\n");
	}
}

my $countuniqSW1 = @uniqSW1;
my $countuniqSW2 = @uniqSW2;

print "Intersection Pipelines File created.\n";
close ($fh);

###################################################################
#Only SW1

my %h;

@h{@uniqSW2}=undef;

my @onlySW1 = grep {not exists $h{$_}} @uniqSW1;

open(my $opp, '>', $onlySW1Perl) or die "cannot open file";

foreach my $op (@onlySW1){
	print($opp "$op\n");
}
 
print "(Only Software 1) Genes File created.\n";
close ($opp);

###################################################################
#Only SW2

my %h1;

@h1{@uniqSW1}=undef;

my @onlySW2 = grep {not exists $h1{$_}} @uniqSW2;


open(my $omp, '>', $onlySW2Perl) or die "cannot open file";

foreach my $om (@onlySW2){
	print($omp "$om\n");
}

print "(Only Software 2) Genes File created.\n";
close ($omp);

######################################################################
#Criação do diagrama de Venn


# Create the Venn::Chart constructor
my $venn_chart = Venn::Chart->new( 550, 550 ) or die("error : $!");
 
# Set a title and a legend for our chart
$venn_chart->set_options( -title => 'Venn diagram - Pipelines Comparison' );
#$venn_chart->set_legends( 'SW1', 'SW2' );
$venn_chart->set_legends( $ARGV[0], $ARGV[2] ); 
# 3 lists for the Venn diagram

 
# Create a diagram with gd object
my $gd_venn = $venn_chart->plot( \@uniqSW1, \@uniqSW2 );
 
# Create a Venn diagram image in png, gif and jpeg format
open my $fh_venn, '>', 'VennChart.png' or die("Unable to create png file\n");
binmode $fh_venn;
print {$fh_venn} $gd_venn->png;
close $fh_venn or die('Unable to close file');
 
# Create an histogram image of Venn diagram (png, gif and jpeg format)
my $gd_histogram = $venn_chart->plot_histogram;
open my $fh_histo, '>', 'VennHistogram.png' or die("Unable to create png file\n");
binmode $fh_histo;
print {$fh_histo} $gd_histogram->png;
close $fh_histo or die('Unable to close file');
 
# Get data list for each intersection or unique region between the 3 lists
my @ref_lists = $venn_chart->get_list_regions();
my $list_number = 1;
foreach my $ref_region ( @ref_lists ) {
 # print "List $list_number : @{ $ref_region }\n";
  $list_number++;
}
print "Venn diagram created.\n";
###################################################################################
#Gererating Reports


my $seqin = Bio::SeqIO->new(-file => "/home/geraldo/Documentos/mestrado/Defesa/GCF_003713225.1_Cara_1.0_genomic.gbff" );

#$seqout= Bio::SeqIO->new( -format => 'Fasta', -file => '>output.fa');

my $reportonlySW2 = "/home/geraldo/Documentos/Fatec/MestradoUTFPR/src/mestradoutfpr/reportonlySW2.txt"; 
my $reportonlySW1 = "/home/geraldo/Documentos/Fatec/MestradoUTFPR/src/mestradoutfpr/reportonlySW1.txt";

open(my $rom, '>', $reportonlySW2) or die "cannot open file";
open(my $rop, '>', $reportonlySW1) or die "cannot open file";



my $reportSW2 = "/home/geraldo/Documentos/Fatec/MestradoUTFPR/src/mestradoutfpr/reportSW2.txt"; 
my $reportSW1 = "/home/geraldo/Documentos/Fatec/MestradoUTFPR/src/mestradoutfpr/reportSW1.txt";

open(my $rm, '>', $reportSW2) or die "cannot open file";
open(my $rp, '>', $reportSW1) or die "cannot open file";



my $reportINTERSEC = "/home/geraldo/Documentos/Fatec/MestradoUTFPR/src/mestradoutfpr/reportintersec.txt";

open(my $ri, '>', $reportINTERSEC) or die "cannot open file";


my $key;
my $content;
my %genomicdata=();

#######################################################################
#Generates all data
my $seq_object = $seqin->next_seq;

 for my $feat_object ($seq_object->get_SeqFeatures) {
    
    for my $tag ($feat_object->get_all_tags) {
    
        for my $value ($feat_object->get_tag_values($tag)) {
			if ($tag eq 'gene'){
				$key = $value;
			}
			if ($tag eq 'note'){
				$content = $value;
			}	
        }
    }
    $genomicdata{$key}=$content;
}


#####################################################################
#Creating report files

print ($rm "\t\t\t\t\tExample of Genes found on Software 2 pipeline\n\n");
	
print ($rp "\t\t\t\t\tExample of Genes found on Software 1 pipeline\n\n");

print ($ri "\t\t\t\t\tExample of Genes found on Both pipelines (Intersection operation)\n\n");

for (keys %genomicdata){
	if(contains($_,\@uniqSW2)){
		print ($rm "Gene: ".$_."\t Obs.:\t".$genomicdata{$_}."\n\n");
	}
}

for (keys %genomicdata){
	if(contains($_,\@uniqSW1)){
		print ($rp "Gene: ".$_."\t Obs.:\t".$genomicdata{$_}."\n\n");
	}
}

for (keys %genomicdata){
	if(contains($_,\@INTERSECTION)){
		print ($ri "Gene: ".$_."\t Obs.:\t".$genomicdata{$_}."\n\n");
	}
}

close $rm;

close $rp;

close $ri;


#####################################################################
#Creating only report files

print ($rom "\t\t\t\t\tExample of Genes only found on Software 2 pipeline\n\n");
	
print ($rop "\t\t\t\t\tExample of Genes only found on Software 1 pipeline\n\n");

for (keys %genomicdata){
	if(contains($_,\@onlySW2)){
		print ($rom "Gene: ".$_."\t Obs.:\t".$genomicdata{$_}."\n\n");
	}
}

for (keys %genomicdata){
	if(contains($_,\@onlySW1)){
		print ($rop "Gene: ".$_."\t Obs.:\t".$genomicdata{$_}."\n\n");
	}
}

close $rom;

close $rop;

print "Report files created\n";

exit 0;
