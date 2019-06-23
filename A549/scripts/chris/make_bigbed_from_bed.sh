#!/bin/bash

bed_list=$(find . -name "*_peaks.bed" | cut -d'/' -f2 | sort)
echo $bed_list

nb_files=$(find . -name "*_peaks.bed" | cut -d'/' -f2 | sort | wc -l)
echo " >>> $nb_files bed files"

hg38_chrom_sizes="/home/chris/Bureau/sb_cofactor_hr/A549/input/hg38.chrom.sizes"

for bed in $bed_list
do
	echo "### $bed"

	echo -e "\t# Sorting..."
	basename="$(cut -d'.' -f1 <<<"$bed")"
	sorted_bed="${basename}.sorted.bed"
	echo -e "\tsort -k1,1 -k2,2n $bed > $sorted_bed"
	# sort -k1,1 -k2,2n $bed > $sorted_bed

	echo -e "\t# Make bigbed..."
	bigbed="${basename}.sorted.bb"
	# bedToBigBed $sorted_bed $hg38_chrom_sizes $bigbed

	echo -e "\t# DONE"
done

# bedToBigBed come from UCSC tools
