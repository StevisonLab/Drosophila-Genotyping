#!/bin/bash

# IMPORTANT: This script is currently hard-coded to Marker 1.
# To make changes for other markers any file with "Marker1" in its name in any command should be simply changed to the appropriate marker number.
# It is probably best to do this with a sed command like sed 's/Marker1/Marker2/g' to be certain all are replaced (this will also replace the string in these instructions)
# Next, the bcftools command below must have genotypes changed appropriately so the logical pattern desired is found

# First bcftools view finds a genotype pattern in the given lines
# The -i option tests if genotypes (GT) are homozygous reference (RR) or homozygous alternate (AA)
# Genotypes are numbered by column number in VCF file
# The -s option is not necessary but can help determine that the given lines match the genotype patterns being asked for
# This is piped to a bcftools query command isolate only the chromosome and position with -f
# Again, the -s option is not necessary but can help validate the command is working as intended, it must include a sample from the -s command in the previous command
# This is piped to awk to print the chromosome number (with "chr" string appended at front to match reference FASTA file)
# awk also will print the position (column $2) with 300 flanking base pairs to make a bed file

bcftools view -s line_42,line_57,line_217,line_357,line_391,line_399,line_437,line_491,line_508,line_810 -i 'GT[8]=="RR" && GT[12]=="RR" && GT[40]=="RR" && GT[79]=="RR" && GT[101]=="RR" && GT[105]=="AA" && GT[111]=="AA" && GT[118]=="AA" && GT[122]=="AA" && GT[174]=="AA"' dgrp2.vcf.gz | bcftools query -s line_42 -f '[%CHROM %POS ]\n' | awk '{OFS="\t"; print "chr"$1, $2-301, $2+300}' > DGRP_Marker1_Candidates.bed

# Then use bedtools to make FASTA files of all sequences defined in bed file

bedtools getfasta -fi dm3.fa -bed DGRP_Marker1_Candidates.bed > DGRP_Marker1_Candidates.fa

# Then define array of sequences excluding sequence names

Sequences=(`grep -v ">" DGRP_Marker1_Candidates.fa`)

# Loop through these sequences and check that the variant falls within the restriction sequence (GAATTC)
# Isolate sequences where this is true
rm DGRP_Marker1_Candidates_2.fa
for Seq in ${Sequences[@]}
do GAATTC=`echo ${Seq:295:11} | grep "GAATTC"`
	if [ -n "${GAATTC}" ]
	then LINE=`grep -n -m 1 "${Seq}" DGRP_Marker1_Candidates.fa | awk -F: '{print $1}'`
	awk -v Line=${LINE} 'NR>Line-2 && NR<Line+1' DGRP_Marker1_Candidates.fa >> DGRP_Marker1_Candidates_2.fa
	fi
done

# Of the above sequences, check that GAATTC is not present multiple times
# Isolate sequences where this is true

COUNTS=(`grep -no "GAATTC" DGRP_Marker1_Candidates_2.fa | uniq -c | awk '{print $1}'`)

ROWS=(`grep -no "GAATTC" DGRP_Marker1_Candidates_2.fa | uniq -c | awk '{print $2}' | awk -F: '{print $1}'`)
rm DGRP_Marker1.fa
for ((LINE=0; LINE<=${#ROWS[@]}-1; LINE++))
do
# index COUNTS and ROWS using [@]:offset:length notation as this is apparently compatible between bash and zsh 
        if [ ${COUNTS[@]:${LINE}:1} -eq 1 ]
	then awk -v ROW=${ROWS[@]:${LINE}:1} 'NR==ROW-1' DGRP_Marker1_Candidates_2.fa >> DGRP_Marker1.fa
	awk -v ROW=${ROWS[@]:${LINE}:1} 'NR==ROW' DGRP_Marker1_Candidates_2.fa >> DGRP_Marker1.fa
	fi
done

# Use seqkit to identify GC content of sequences and rank them by highest gc content
# This requires seqkit to be downloaded on computer and isn't strictly necessary as GC concents of primer regions can be validated separately

seqkit fx2tab --gc DGRP_Marker1.fa | awk '{print $1, $3}' | sort -r -k2,2 > DGRP_Marker1_List.txt
