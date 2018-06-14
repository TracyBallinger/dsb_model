#!/bin/bash 

# Download data from "High-throughput sequencing of DNA G-quadruplex structures in the human genome" by Chambers et al. Nature Biotechnology 2015
LC_COLLATE=C
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874_Na_PDS_12_minus.bedGraph.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874_Na_PDS_12_plus.bedGraph.gz
# Add the plus and minus strand signals together.
bedtools unionbedg -i GSE63874_Na_PDS_12_minus.bedGraph.gz GSE63874_Na_PDS_12_plus.bedGraph.gz > testout
awk '{print $1"\t"$2"\t"$3"\t"$4+$5}' testout > GSE63874_Na_PDS_12_both.bedGraph

chromsizes=hg19_chr.lengths
bedGraphToBigWig GSE63874_Na_PDS_12_both.bedGraph $chromsizes GSE63874_Na_PDS_12_both.bigWig

