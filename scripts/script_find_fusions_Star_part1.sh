#!/bin/bash

## Script to create the reference file for each translocation

## input files:
# file with the paths for Ensembl annotation and genome files: "path_file.txt"
# file with the list of translocations to be searched: "translocation_list.txt"
# 
#"path_file.txt":
#1 row: path to annotation file
#2 row: path to genome file
# 
#"translocation_list.txt":
#1 field: gene name in 5'
#2 field: gene name in 3'
#3 field: exon number of the gene in 5’ where the direct breakpoint occurs (optional)
#4 field: exon number of the gene in 3’ where the direct breakpoint occurs (optional)
#5 field: exon number of the gene in 3’ where the reciprocal breakpoint occurs (optional)
#6 field: exon number of the gene in 5’ where the reciprocal breakpoint occurs (optional)
#7 field: Ensembl id for the gene in 5’
#8 field: Ensembl id for the gene in 3’

## command: /scripts/script_find_fusions_Star_part1.sh gene_name1 gene_name2



# input
GENE1=$1
GENE2=$2

file_annotation=$( cat /data/path_file.txt | sed -n '1 p' )	# "../../annotation/Homo_sapiens.GRCh38.97.gtf"
file_genome=$( cat /data/path_file.txt | sed -n '2 p' )	# "../../annotation/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
file_translocation="/data/translocation_list.txt"



# 1) Retrieve exon regions of fusion genes

# gene1

START1=$( cat $file_annotation | awk '$3=="gene" {print}' | grep -w "\"$GENE1\"" | cut -f4 )
END1=$( cat $file_annotation | awk '$3=="gene" {print}' | grep -w "\"$GENE1\"" | cut -f5 )
CHR1=$( cat $file_annotation | awk '$3=="gene" {print}' | grep -w "\"$GENE1\"" | cut -f1 )
STRAND1=$( cat $file_annotation | awk '$3=="gene" {print}' | grep -w "\"$GENE1\"" | cut -f7 )

echo $GENE1 $START1 $END1 $CHR1 $STRAND1 | sed -e 's/ /\t/g' > info_$GENE1.txt

TRANSCRIPT1=$( cat $file_translocation | awk '$1=="'$GENE1'" && $2=="'$GENE2'" {print $7}' )

cat $file_annotation | grep -w "\"$GENE1\"" | grep -w "\"$TRANSCRIPT1\"" | grep -w "exon" | nl | awk '{print $2, $5-1, $6, "exon_"$1, "'$GENE1'", $8}' | sed -e 's/ /\t/g' > exons_$GENE1.bed


# gene2

START2=$( cat $file_annotation | awk '$3=="gene" {print}' | grep -w "\"$GENE2\"" | cut -f4)
END2=$( cat $file_annotation | awk '$3=="gene" {print}' | grep -w "\"$GENE2\"" | cut -f5)
CHR2=$( cat $file_annotation | awk '$3=="gene" {print}' | grep -w "\"$GENE2\"" | cut -f1)
STRAND2=$( cat $file_annotation | awk '$3=="gene" {print}' | grep -w "\"$GENE2\"" | cut -f7)

echo $GENE2 $START2 $END2 $CHR2 $STRAND2 | sed -e 's/ /\t/g' > info_$GENE2.txt

TRANSCRIPT2=$( cat $file_translocation | awk '$1=="'$GENE1'" && $2=="'$GENE2'" {print $8}' )

cat $file_annotation | grep -w "\"$GENE2\"" | grep -w "\"$TRANSCRIPT2\"" | grep -w "exon" | nl | awk '{print $2, $5-1, $6, "exon_"$1, "'$GENE2'", $8}' | sed -e 's/ /\t/g' > exons_$GENE2.bed





# 2) Retrieve exon sequences of the fusion genes

bedtools getfasta -fi $file_genome -bed exons_$GENE1.bed -s -name -tab -fo exons_$GENE1.txt

bedtools getfasta -fi $file_genome -bed exons_$GENE2.bed -s -name -tab -fo exons_$GENE2.txt





# 3) Find all the possible fusion combinations between the exons of the two fusion genes

mv exons_$GENE1.txt exons_gene1.txt
mv exons_$GENE2.txt exons_gene2.txt

R CMD BATCH /scripts/script_crossjoin_per_fusions.r

mv exons_gene1.txt exons_$GENE1.txt
mv exons_gene2.txt exons_$GENE2.txt

mv exons_forward.txt exons_$GENE1-$GENE2.txt
mv exons_reverse.txt exons_$GENE2-$GENE1.txt





# 4) Create the reference with the breakpoint position

# direct translocation:

# retrieve the breakpoint position
cat exons_$GENE1-$GENE2.txt | awk '{print $1, "-", $3}' | sed -e 's/ //g' > a1.txt
cat exons_$GENE1-$GENE2.txt | cut -f2 | awk '{print length}' | awk '{print $1-1, $1+1}' | sed -e 's/ /\t/' > b1.txt
paste a1.txt b1.txt > reference_$GENE1-$GENE2.bed

# create the reference
cat exons_$GENE1-$GENE2.txt | awk '{print "'$GENE1'", "-", "'$GENE2'", "|", $1, "-", $3, $2, $4}' | sed -e 's/ /\t/7; s/ //g' > exons_for_align_$GENE1-$GENE2.txt
cat exons_for_align_$GENE1-$GENE2.txt | sed -e 's/^/>/; s/\t/\n/' > exons_for_align_$GENE1-$GENE2.fa

# remove unuseless files
rm a1.txt
rm b1.txt

# add each single reference sequence to the all-encompassing reference sequence file
cat exons_for_align_$GENE1-$GENE2.fa >> /data/reference/reference_sequences.fa



# reciprocal translocation:

# retrieve the breakpoint position
cat exons_$GENE2-$GENE1.txt | awk '{print $1, "-", $3}' | sed -e 's/ //g' > a2.txt
cat exons_$GENE2-$GENE1.txt | cut -f2 | awk '{print length}' | awk '{print $1-1, $1+1}' | sed -e 's/ /\t/' > b2.txt
paste a2.txt b2.txt > reference_$GENE2-$GENE1.bed

# create the reference
cat exons_$GENE2-$GENE1.txt | awk '{print "'$GENE2'", "-", "'$GENE1'", "|", $1, "-", $3, $2, $4}' | sed -e 's/ /\t/7; s/ //g' > exons_for_align_$GENE2-$GENE1.txt
cat exons_for_align_$GENE2-$GENE1.txt | sed -e 's/^/>/; s/\t/\n/' > exons_for_align_$GENE2-$GENE1.fa

# remove unuseless files
rm a2.txt
rm b2.txt

# add each single reference sequence to the all-encompassing reference sequence file
cat exons_for_align_$GENE2-$GENE1.fa >> /data/reference/reference_sequences.fa



# remove unuseless files
rm script_crossjoin_per_fusions.r.Rout
rm .RData

