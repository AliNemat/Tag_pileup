# Tag_pileup

A bioinformatics tool for tag pileup



## The python code takes four inputs:

1- A gff file *

2- Sam file **

3-Desired reference point for tagpile up ***

4-The length of first-pair read in pair end sequencing ****


*Each row of gff file provides chromosome, first bp genome location of the specific motif, last bp genome location of the specific motif, motif number (1, 2 or 3 if three motifs are pogrammed to find), strand direction (positive or negative) 

** Sam file is typically deduplicated mapped read

Information this code extract from SAM file are: 

a) Number of tags in 5' end of first pair read (Second pair and rest of first pair read won't be used in chipexo assay)

b) Only reads which are maaped with map quality > 20 (99% accuracy) will be considered.

c) Tags belong to positive or negative strand of DNA



*** Example: If motif is  CCGTTACGATCGA and we wish the first A be our reference point, we should enter
   6

****  Example: If sequencer is set to read 40 bp for the first end pair and 36 for the second pair, 40 needs to be given as input.

Why do we need this information?

In the positive strand of DNA, the left-most position of the first pair read, reported in SAM file, is identical to 5' end.
However, in the negative strand of DNA, left-most position of the first pair read plus length of first pair read gives 5' end of the read.


## R code
The R code is a simple postprocessor to get statistics, plots and compare tagpile up results   


Tagpile usage in a typical pipeline: 

A peakcaller tool (for example: Genetrak/CW-PAir, Chexmix) finds specific loctions in genome as preaks (events)
This locations will be expanded (for example -30bp, 30bp) and filtered (For example by RepeatMasker) and then feeds into MEME to finds the motif candidates.

Output of MEME which are specific letters and weights are given to FIMO along with expanded for example -30, 30) but not filtered genome locations to search for motifs predicted by MEME.

Ouput put of FIMO in gff format is given as input to Tagpileup along with the SAM file.
