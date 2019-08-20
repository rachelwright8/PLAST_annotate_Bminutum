# PLAST_annotate_Bminutum

BLASTx takes TOO LONG!

I'm giving PLAST a shot instead.

https://plast.inria.fr/plast-algorithm/

Other users have found PLAST to be way faster and nearly as sensitive as BLAST.

https://2-bitbio.com/2017/07/blastx-is-too-slow-heres-some.html

This is my walkthru for annotating the Breviolum minutum transcriptome.

I was able to run the annotations on my personal computer (iMac, 3 GHz Intel Core i5, 16 GB memory) in just a couple of hours.

I compared my annotations to BLAST annotations made by Dr. John Parkinson. See the 'test' directory for an R script to compare the annotations.

In summary: 
PLAST runs ways faster than BLAST.
I get nearly as many genes and the genes are identical to BLAST.
The GO terms found differed more, but this may be due to updates to the uniprot database since the BLAST annotation was peformed (several years prior).
