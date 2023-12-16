#!/bin/bash

## Prepare workspace

# activate the ```starfish``` conda environment:
# conda activate starfish

# copy the ```starfish``` example data directory into a test directory:
cp -r examples/ test/
cd test/

# create ome2\*.txt files detailing the absolute path to each genome's gff3 and assembly:
realpath assembly/* | perl -pe 's/^(.+?([^\/]+?).fasta)$/\2\t\1/' > ome2assembly.txt
realpath gff3/* | perl -pe 's/^(.+?([^\/]+?).final.gff3)$/\2\t\1/' > ome2gff.txt

# concatenate all gff3 files into a single file (a useful shortcut for some analyses):
cat gff3/*.gff3 > macpha6.gff3

# concatenate all assembly files and make a blastn database:
mkdir blastdb
cut -f2 ome2assembly.txt | xargs cat > blastdb/macpha6.assemblies.fna
makeblastdb -in blastdb/macpha6.assemblies.fna -out blastdb/macpha6.assemblies -parse_seqids -dbtype nucl

# calculate %GC content across all genomes (useful for visualizing elements later):
../aux/seq-gc.sh -Nbw 1000 blastdb/macpha6.assemblies.fna > macpha6.assemblies.gcContent_w1000.bed
rm blastdb/macpha6.assemblies.fna

# parse the provided eggnog mapper annotations (NB the format of the output file has changed in more recent emapper versions):
cut -f1,12  ann/*emapper.annotations | grep -v  '#' | grep -v -P '\t-' | perl -pe 's/\t/\tEMAP\t/' | grep -vP '\tNA' > ann/macph6.gene2emap.txt

# retrieve the narrowest eggnog ortholog group per sequence and convert to mcl format:
cut -f1,10 ann/*emapper.annotations | grep -v '#' | perl -pe 's/^([^\s]+?)\t([^\|]+).+$/\1\t\2/' > ann/macph6.gene2og.txt

# convert to .mcl format:
../aux/geneOG2mclFormat.pl -i ann/macph6.gene2og.txt -o ann/

## Gene finder module

# We begin by *de novo* annotating all tyrosine recombinases (tyrs/YRs) in the provided assemblies. In practice, we can *de novo* annotate any gene we want, as long as we have an HMM file of a predicted domain within that gene and a multifasta of amino acid sequences of that gene (the more predicted sequences the better).

# first, create a dedicated directory for good housekeeping:
mkdir geneFinder

# *de novo* annotate tyrs with the provided YR HMM and amino acid queries (~10min):
starfish annotate -T 2 -x macpha6_tyr -a ome2assembly.txt -g ome2gff.txt -p ../db/YRsuperfams.p1-512.hmm -P ../db/YRsuperfamRefs.faa -i tyr -o geneFinder/

# you should observe the following printed output:
# found 1 new tyr genes and 9 tyr genes that overlap with 11 existing genes

# now consolidate the newly predicted gene coordinates with the existing gff3:
starfish consolidate -o ./ -g macpha6.gff3 -G geneFinder/macpha6_tyr.filt_intersect.gff

# create a .txt file with the path to the new consolidated gff file:
realpath macpha6_tyr.filt_intersect.consolidated.gff | perl -pe 's/^/macpha6\t/' > ome2consolidatedGFF.txt

# organize tyrs into mutually exclusive neighbourhoods separated by at least 10kb to avoid adjacent tyrs messing up subsequent analyses:
starfish sketch -m 10000 -q geneFinder/macpha6_tyr.filt_intersect.ids -g ome2consolidatedGFF.txt -i s -x macpha6 -o geneFinder/

# you should observe the following printed output:
# found 10 neighborhoods containing input query genes

# neighbourhoods will often contain intervening genes located between genes of interest so pull out the coordinates of candidate captains only:
grep -P '\ttyr\t' geneFinder/macpha6.bed > geneFinder/macpha6.tyr.bed 

#### Key output files

# * the \*.filt_intersect.\* files contain all newly predicted genes and amino acid sequences. Newly predicted genes that overlap with an existing gene will keep their original sequenceID (or if overlapping with multiple existing genes, will be assigned a new sequenceID consisting of the concatenated existing gene IDs) but will be assigned the newly predicted amino acid sequence. Newly predicted genes that don't overlap with an existing gene will be assigned a new sequenceID. 
# * the *.bed file contains the coordinates of all genes of interest organized into neighbourhoods

## Element finder module

# Now we can annotate mobile elements containing the candidate captain genes. In order to be found, elements must have the basic architecture of a fungal *Starship* or bacterial integrative and conjugative element: a captain gene with zero or more cargo genes downstream of its 3' end.

# create a dedicated directory for good housekeeping:
mkdir elementFinder

# search for insertions containing at least one predicted tyr (takes ~1min):
starfish insert -T 2 -a ome2assembly.txt -d blastdb/macpha6.assemblies -b geneFinder/macpha6.tyr.bed -i tyr -x macpha6 -o elementFinder/

# you should observe the following printed output:
# found element boundaries and insertion sites for 7 tyr captains out of 10 input captains

# search for flanking repeats around predicted element boundaries:
starfish flank -a ome2assembly.txt -b elementFinder/macpha6.insert.bed -x macpha6 -o elementFinder/

# you should observe the following printed output:
# found 2 captains with DR boundaries and 5 captains with DR-TIR boundaries out of 7 input captains with candidate DRs

# summarize the element metadata, identify overlaps, name sites and identify all captain and cargo genes:
starfish summarize -a ome2assembly.txt -b elementFinder/macpha6.flank.bed -x macpha6 -o elementFinder/ -S elementFinder/macpha6.insert.stats -f elementFinder/macpha6.flank.singleDR.stats -g ome2consolidatedGFF.txt -A ann/macph6.gene2emap.txt -t geneFinder/macpha6_tyr.filt_intersect.ids 

# it is strongly recommended to look at an alignment of each element against its 'best' insertion site to manually filter out false positives. Use circos to visualize nucmer alignments (takes ~2min):
mkdir pairViz
starfish pair-viz -m all -t empty -T 2 -A nucmer -a ome2assembly.txt -b elementFinder/macpha6.elements.bed -f elementFinder/macpha6.flank.singleDR.stats -S elementFinder/macpha6.elements.named.stats -o pairViz/

#### Key output files

# * \*.insert.bed contains coordinates of all predicted element boundaries based on all candidate insertions. 
# * \*.insert.stats contains useful metadata on candidate insertions.
# * \*.flank.bed contains updated coordinates of predicted element boundaries based on where flanking repeats were found. Alternatively, if no flanking repeats were found, it contains the coordinates of the boundaries that give the longest element. 
# * \*.flank.singleDR.stats contains useful metadata on recovered flanking repeats.
# * \*.elements.bed contains all captain, boundary, and gene features of predicted elements. 
# * \*.elements.feat contains element metadata. 
# * \*.elements.fna contains element sequences. 
# * \*.elements.named.stats is an updated version of the insert.stats file that contains named insertion sites.

## Region finder module

# The final step of a ```starfish``` analysis involves situating elements and insertion sites into homologous genomic regions so we know whether insertions are shared across individuals or not. We leverage two types of orthologous relationships: between elements and between all genes across all genomes.

# You can group elements together based solely on the orthology of their captains. However, since captains from the same ortholog family often have different cargo haplotypes, I recommend using a combination of captain families and cargo haplotypes to group elements by "family:haplotype" combinations. 

# We use ortholog information from the eggnog mapper analysis to group all genes into ortholog groups since this is just a quick analysis. I would otherwise recommend running a more comprehensive analysis, like Orthofinder.

# create a dedicated directory for good housekeeping:
mkdir regionFinder

# group all tyrs into families using ```mmseqs2 easy-clust``` with a very permissive 50% percent ID/ 25% coverage threshold (families with only a single member will automatically be assigned the prefix 'sng'):
mmseqs easy-cluster geneFinder/macpha6_tyr.filt_intersect.fas regionFinder/macpha6_tyr regionFinder/ --threads 2 --min-seq-id 0.5 -c 0.25 --alignment-mode 3 --cov-mode 0 --cluster-reassign
../aux/mmseqs2mclFormat.pl -i regionFinder/macpha6_tyr_cluster.tsv -g fam -o regionFinder/

# use ```sourmash``` and ```mcl``` to group all elements into haplotypes based on pairwise k-mer similarities across entire elements:
starfish sim -m element -t nucl -b elementFinder/macpha6.elements.bed -x macpha6 -o regionFinder/ -a ome2assembly.txt
starfish group -m mcl -s regionFinder/macpha6.element.nucl.sim -i hap -o regionFinder/ -t 0.05

# replace captainIDs with elementIDs in the captain groups file:
grep -P '\tcap\t' elementFinder/macpha6.elements.bed | cut -f4,7 > regionFinder/macpha6.cap2ship.txt
../aux/searchReplace.pl -i regionFinder/macpha6_tyr_cluster.mcl -r regionFinder/macpha6.cap2ship.txt > regionFinder/macpha6.element_cluster.mcl

# merge captain family with element haplotype info:
../aux/mergeGroupfiles.pl -t regionFinder/macpha6.element_cluster.mcl -q regionFinder/macpha6.element.nucl.I1.5.mcl > regionFinder/macpha6.element.fam-hap.mcl

# create a file with tyrs that are not found in any elements (will let us assign them to fragmented haplotypes in the ```dereplicate``` analysis):
grep -f <(comm -23 <(cut -f1 geneFinder/macpha6_tyr.filt_intersect.ids | sort) <(grep -P '\tcap\t|\ttyr\t' elementFinder/macpha6.elements.bed | cut -f4| sort)) geneFinder/macpha6.tyr.bed > regionFinder/unaffiliated_tyrs.bed

# you can increase confidence in region homology by only looking at gene ortholog groups with low copy numbers missing from few genomes:
../aux/filterOG.pl -O ann/macph6.gene2og.mcl -a 1 -c 5 -o ann/

# It is more useful to play around with copy number thresholds than with genome absence thresholds because ```dereplicate``` will automatically filter OGs to retain those that have more taxonomic information i.e., are present in a greater number of individuals.

# now, dereplicate your data:
starfish dereplicate -e regionFinder/macpha6.element.fam-hap.mcl -t regionFinder/unaffiliated_tyrs.bed -F elementFinder/macpha6.elements.feat -S elementFinder/macpha6.elements.named.stats -O ann/macph6.gene2og.a1.c5.txt -g ome2gff.txt -x macpha6 -o regionFinder/ --flanking 3

# I would normally recommend going with the default ```--flanking 6``` but because the *Defiant* insertion is in a gene-sparse region, it can only be recovered with ```---flanking 3```

# a lot of information will be printed to STDOUT, but somewhere in there you should see:
# found 5 regions with at least 1 cross-referenced element-insertion site pair

# Element haplotypes consist of predicted mobile element sequences. Empty haplotypes consist of a contiguous sequence formed by the flanking regions of an element. Fragmented haplotypes consist of a non-empty sequence flanked by the flanking regions of an element but missing a predicted element. Coordinates of predicted insertion sites from the elementFinder module are cross referenced with empty and fragmented haplotypes and are considered to be 'verified' if their coordinates overlap. 

# as before, it is strongly recommended to look at haplotype alignments within each region to manually filter out false positives. Use gggenomes to visualize nucmer alignments (takes ~5min):
mkdir locusViz
starfish locus-viz -T 2 -m region -a ome2assembly.txt -b elementFinder/macpha6.elements.bed -x macpha6 -o locusViz/ -A nucmer -r regionFinder/macpha6.fog3.d600000.m1.regions.txt -d regionFinder/macpha6.fog3.d600000.m1.dereplicated.txt -j regionFinder/macpha6.fog3.d600000.m1.haplotype_jaccard.sim  -g ome2consolidatedGFF.txt --tags geneFinder/macpha6_tyr.filt_intersect.ids --gc macpha6.assemblies.gcContent_w1000.bed