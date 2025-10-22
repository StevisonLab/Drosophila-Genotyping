#!/bin/bash

# IMPORTANT: This script is currently hard-coded to Marker 1.
# To make changes for other markers any file with "Marker1" in its name in any command should be simply changed to the appropriate marker number.
# It is probably best to do this with a sed command like sed 's/Marker1/Marker2/g' to be certain all are replaced (this will also replace the string in these instructions)
# Next, the bcftools command below must have genotypes changed appropriately so the logical pattern desired is found

# First bcftools view finds a genotype pattern in the given lines
# The -i option tests if genotypes (GT) are homozygous reference (RR) or homozygous alternate (AA)
# Genotypes are numbered by column number in VCF file
# This is piped to a bcftools query command to isolate only the chromosome and position with -f
# The -s option is not necessary but can help validate the command is working as intended, it must include a sample from the -s command in the previous command
# This is piped to awk to print the chromosome number (with "chr" string appended at front to match reference FASTA file)
# awk also will print the position (column $2) with 300 flanking base pairs to make a bed file
# The -i command here also includes conditions for the Genotyping Quality (phred score) and read depth (sometimes called coverage) to be above certain thresholds. These can be changed if needed.

bcftools view -i 'MIN(FORMAT/GQ)>20 && MIN(FORMAT/DP)>10 && GT[0]=="AA" && GT[1]=="RR" && GT[2]=="AA" && GT[3]=="AA" && GT[4]="RR"' Filters.vcf.gz | bcftools query -s Flagstaff16 -f '[%CHROM %POS ]\n' | awk '{OFS="\t"; print $1, $2-301, $2+300}' > DPSE_Marker1_Candidates.bed

# Then use bedtools to make FASTA files of all sequences defines in bed file

bedtools getfasta -fi GCF_009870125.1_UCI_Dpse_MV25_genomic.fna -bed DPSE_Marker1_Candidates.bed > DPSE_Marker1_Candidates.fa

# Then define array of sequences excluding sequence names

Sequences=(`grep -v ">" DPSE_Marker1_Candidates.fa`)

# Loop through these sequences and check that the variant falls within the restriction sequence (GAATTC)
# Isolate sequences where this is true

for Seq in ${Sequences[@]}
do GAATTC=`echo ${Seq:295:11} | grep "GAATTC"`
	if [ -n "${GAATTC}" ]
	then LINE=`grep -n -m 1 "${Seq}" DPSE_Marker1_Candidates.fa | awk -F: '{print $1}'`
	awk -v Line=${LINE} 'NR>Line-2 && NR<Line+1' DPSE_Marker1_Candidates.fa >> DPSE_Marker1_Candidates_2.fa
	fi
done

# Of the above sequences, check that GAATTC is not present multiple times
# Isolate sequences where this is true

COUNTS=(`grep -no "GAATTC" DPSE_Marker1_Candidates_2.fa | uniq -c | awk '{print $1}'`)

ROWS=(`grep -no "GAATTC" DPSE_Marker1_Candidates_2.fa | uniq -c | awk '{print $2}' | awk -F: '{print $1}'`)

for ((LINE=0; LINE<=${#ROWS[@]}-1; LINE++))
do
	if [ ${COUNTS[@]:${LINE}:1} -eq 1 ]
	then awk -v ROW=${ROWS[@]:${LINE}:1} 'NR==ROW-1' DPSE_Marker1_Candidates_2.fa >> DPSE_Marker1.fa
	awk -v ROW=${ROWS[@]:${LINE}:1} 'NR==ROW' DPSE_Marker1_Candidates_2.fa >> DPSE_Marker1.fa
	fi
done

# Use seqkit to identify GC content of sequences and rank them by highest gc content
# This requires seqkit to be downloaded on computer and isn't strictly necessary as GC concents of primer regions can be validated separately

seqkit fx2tab --gc DPSE_Marker1.fa | awk '{print $1, $3}' | sort -r -k2,2 > DPSE_Marker1_List.txt
