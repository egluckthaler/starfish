CNEFinder: Finding Conserved Non-coding Elements in Genomes
===

NOTE: this version of CNEFinder has been modified to analyze a minimum repeat length of 4 so that it can be used by starfish, the Starship Finder Shell. More information available at https://github.com/egluckthaler/starfish

<b>Description</b>: Given two genomes r and q, and a reference and query gene (or coordinate values), CNEFinder finds CNE matches between either the positions of the provided genes or coordinate values of the respective chromosomes of r and q.

<b>Installation</b>: To compile CNEFinder, please follow the instructions given in file INSTALL.
```
 Usage: CNEFinder <options>
 Standard (Mandatory):
  -r, --ref-genome-file		<str>		FASTA reference genome filename.
  -q, --query-genome-file	<str>		FASTA query genome filename.
  -e, --exons-ref-file		<str>		TSV exon coordinates for reference genome filename.
  -f, --exons-query-file	<str>		TSV exon coordinates for query genome filename.
  -l, --min-seq-length		<int>		Minimum length of CNE.
  -t, --sim-threshold		<dbl>		Threshold of similarity between sequences (0-1].
  -o, --output-file		<str>		Output filename with CNEs identified in BED format.

  Either 1. or 2.
    1.Search using gene name:
    -g, --ref-gene-file		<str>		TSV filename containing gene data for reference genome.
    -n, --ref-gene-name		<str>		Name of gene in reference genome in which CNEs will be identified.
    -j, --query-gene-file	<str>		TSV filename containing gene data for query genome.
    -m, --query-gene-name	<str>		Name of gene in query genome in which CNEs will be identified.

    2.Search using index position:
    -y, --ref-chrom		<str>		Chromosome of reference genome.
    -z, --query-chrom		<str>		Chromosome of query genome.
    -a, --ref-start		<int>		Start CNE search from this position of reference sequence.
    -b, --ref-end		<int>		End CNE search at this position of reference sequence.
    -c, --query-start		<int>		Start CNE search from this position of query sequence.
    -d, --query-end		<int>		End CNE search at this position of query sequence.

  Optional:
  -Q, --mem-length		<int>		Minimum length of maximal exact matches. Default:18.
  -M, --merged-length		<dbl>		Minimum length (in terms of CNE length) of merged matches to be extended. Default:0.5.
  -s, --ext-threshold		<dbl>		Threshold to further extend similarity threshold by. Default:0.05.
  -u, --max-seq-length		<int>		Set a maximum length for the CNE. Default:2000.
  -p, --repeat-regions		<int>		Choose 1 to filter repetitive regions of genomes or 0 otherwise. Default:1.	
  -v, --rev-complement		<int>		Choose 1 to compute CNEs for reverse complement or 0 otherwise. Default:0.
  -x, --remove-overlaps		<int>		Choose 1 to remove overlapping CNEs or 0 otherwise. Default:1.

 Number of threads:
  -T, --threads			<int>		Number of threads to use. Default:1. 
```

<b>See https://github.com/lorrainea/CNEFinder/wiki for more help.</b>

<b>Citation</b>:

```
L. A. K. Ayad, S. P. Pissis, D. Polychronopoulos: 
CNEFinder: finding conserved non-coding elements in genomes. 
Bioinformatics: 34(17): i743-i747 (2018) 
```
<b>License</b>: GNU GPLv3 License; Copyright (C) 2017 Lorraine A. K. Ayad, Solon P. Pissis, Dimitris Polychronopoulos

