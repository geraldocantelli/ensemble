# ensemble
Ensemble Solution for Comparison of Genomic Annotation Softwares

Once you want to compare two genomic annotation softwares and already executed both pipelines, this tool provides a way to perform it.

All implemented in Perl language with BioPerl support.

Before execute the scripts, run:

   cpan install Venn::Chart
  
   cpan install Bio::SeqIO

Use the scripts called inventorySW1.pl and inventorySW2.pl to make the inventory of the genes resulted of each pipeline, as follows:
If the pipeline process ended in a GFF3 file, use:

   perl ./inventorySW1.pl ./nameofPipelineresultantFile.gff3 ./recordDestSW1.txt

Observe that there are just two parameters: the first is the file that resulted from the pipeline and the second is the name of the file generated by the inventory.

If the pipeline process ended in a GFF file, use:

   perl ./inventorySW2.pl ./nameofPipelineresultantFile.gff ./recordDestSW2.txt

Each script is adapted for the characteristics of each annotation format and extract the gene names, generating an txt file.

Now you have the inventory gene file of each sofware, rename them properly in order to use as entry of ensemble script:
   
   perl ./ensemble.pl nameofSoftware1 ./recordDestSW1.txt nameofSofware2 ./recordDestSW2.txt
   
The ensemble script generates a Venn diagram with labels given by the names of software of the parameters and also five reports in file format:

   onlySW1.txt:     information about the genes present in the pipeline of software 1
  
   onlySW2.txt:     information about the genes present in the pipeline of software 2
  
   reportonlySW1.txt: information about the genes present exclusively in the pipeline of software 1
   
   reportonlySW2.txt: information about the genes present exclusively in the pipeline of software 2
   
   intersec.txt:      information about the genes present in the intersection set of the pipelines

The two diagrams generated are called:

   VennChart.png
   
   VennHistogram.png
