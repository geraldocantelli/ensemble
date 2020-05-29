use warnings;
use Carp;
use strict;
 
use Venn::Chart;
use Array::Contains;
use Array::Utils qw(:all);
use List::Uniq ':all';
use Bio::Seq;
use Bio::SeqIO;


#!/usr/bin/perl
my $gravaDestinoSW1 = $ARGV[1];	
my $gravaDestinoSW2 = $ARGV[3];
my $INTERSEC = "./INTERSEC.txt";

my $fileonlySW1 = "./onlySW1.txt";
my $fileonlySW2 = "./onlySW2.txt";

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

my @uniqSW1= uniq(@SW1);
my @uniqSW2= uniq(@SW2);

###################################################################
#Only SW1

my %h;

@h{@uniqSW2}=undef;

my @onlySW1 = grep {not exists $h{$_}} @uniqSW1;

open(my $opp, '>', $fileonlySW1) or die "cannot open file";

foreach my $op (@onlySW1){
	print($opp "$op\n");
}
 
print "(Only ".$ARGV[0].") Genes File created.\n";
close ($opp);

###################################################################
#Only SW2

my %h1;

@h1{@uniqSW1}=undef;

my @onlySW2 = grep {not exists $h1{$_}} @uniqSW2;


open(my $omp, '>', $fileonlySW2) or die "cannot open file";

foreach my $om (@onlySW2){
	print($omp "$om\n");
}

print "(Only ".$ARGV[2].") Genes File created.\n";
close ($omp);


################################################################################
#Creating Intersection set



my @INTERSECTION = intersect(@onlySW1, @onlySW2);

foreach my $p (@INTERSECTION){
	print($fh "$p\n");

}


print "Intersection File created.\n";
close ($fh);

######################################################################
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

my $reportonlySW2 = "./reportonlySW2.txt"; 
my $reportonlySW1 = "./reportonlySW1.txt";

open(my $rom, '>', $reportonlySW2) or die "cannot open file";
open(my $rop, '>', $reportonlySW1) or die "cannot open file";



my $reportINTERSEC = "./reportintersec.txt";

open(my $ri, '>', $reportINTERSEC) or die "cannot open file";


my $key;
my $content;
my %genomicdata=();
my $gene="";
my $note="";

#######################################################################
#Generates all data



my $seq_object = $seqin->next_seq;

for my $feat_object ($seq_object->get_SeqFeatures) {
    if ($feat_object->primary_tag eq "gene") {
        
        if ($feat_object->has_tag("gene")) {
            for my $val ($feat_object->get_tag_values("gene")) {
                #print $val." ";
		$gene=$val;
	    }
        }
        if ($feat_object->has_tag("note")) {
            for my $val ($feat_object->get_tag_values("note")) {
                #print $val."\n";
		$note=$val;
            }
        }
	$genomicdata{$gene}=$note;
    }
}

#####################################################################
#Creating report files

print ($rom "\t\t\t\t\tExample of Genes found on ".$ARGV[2]." pipeline\n\n");
	
print ($rop "\t\t\t\t\tExample of Genes found on ".$ARGV[0]." pipeline\n\n");

print ($ri "\t\t\t\t\tExample of Genes found on Both pipelines (Intersection operation)\n\n");

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

for (keys %genomicdata){
	if(contains($_,\@INTERSECTION)){
		print ($ri "Gene: ".$_."\t Obs.:\t".$genomicdata{$_}."\n\n");
	}
}

close $rom;

close $rop;

close $ri;

print "Report files created\n";

exit 0;
