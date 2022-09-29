#!/bin/bash

## Script to detect fusion transcripts for each sample

## input files:
# file with the parameters used by Star embedded in CircFusion: "params_Star.txt"
# file with the parameters used by CircFusion to define fusion transcripts: "params_fusions.txt"
# 
# "params_Star.txt"
#1 row: number of threads for the analysis (default 8)
#2 row: uncompression command; the default is “zcat” for gzip compressed file (.gz)
#3 row: maximum number of multiple alignments allowed (default 10)
#4 row: minimum number of matches required for the alignment (default 5)
#5 row: maximum number of mismatches allowed in the alignment (default 4)
#6 row: maximum number of mismatches allowed per pair relative to read length (default 0.01)
#7 row: maximum number of mismatches allowed for stitching of the splice junctions (default -1 -1 -1 -1)
# 
# "params_fusions.txt"
#1 row: maximum percentage of mismatches allowed at the size of shorter overlap over the fusion junction (default 17)
#2 row: minimum number of nucleotides overlapping the fusion junction (default 5)
#3 row: minimum number of nucleotides mapping in the alignment (default 50)
#4 row: library strand; allowed values are "+", "-" or "UN", for stranded, reverse stranded or unstranded RNASeq, respectively (default “+”)

## command: /scripts/script_find_fusions_Star_part2.sh $M $END $SAMPLE $SAMPLE_LIST



M=$1
END=$2
SAMPLE=$3
SAMPLE_LIST=$4


# from "params_Star.txt":
THR=$( cat /data/params_Star.txt | sed -n '1 p' )					# default: 8
UNCOMPRESSION_COMMAND=$( cat /data/params_Star.txt | sed -n '2 p' )			# default: "zcat"
OUTFILTERMULTIMAPNMAX=$( cat /data/params_Star.txt | sed -n '3 p' )			# max num multiple alignments; default: 10
OUTFILTERMATCHNMIN=$( cat /data/params_Star.txt | sed -n '4 p' )			# min num match; default: 5
OUTFILTERMISMATCHNMAX=$( cat /data/params_Star.txt | sed -n '5 p' )			# max num mismatch; default: 4
OUTFILTERMISMATCHNOVERREADLMAX=$( cat /data/params_Star.txt | sed -n '6 p' )		# max num mismatches per pair relative to read length; default: 0.01 -> 2m per 2*100nt, 3m per 2*150nt
ALIGNSJSTITCHMISMATCHNMAX=$( cat /data/params_Star.txt | sed -n '7 p' )		# maximum number of mismatches for stitching of the splice junctions; default: -1 -1 -1 -1

#default values
if [ -z $THR ]
then
	THR=8
fi

if [ -z $UNCOMPRESSION_COMMAND ]
then
	UNCOMPRESSION_COMMAND="\"zcat\""
fi

if [ -z $OUTFILTERMULTIMAPNMAX ]
then
	OUTFILTERMULTIMAPNMAX=10
fi

if [ -z $OUTFILTERMATCHNMIN ]
then
	OUTFILTERMATCHNMIN=5
fi

if [ -z $OUTFILTERMISMATCHNMAX ]
then
	OUTFILTERMISMATCHNMAX=4
fi

if [ -z $OUTFILTERMISMATCHNOVERREADLMAX ]
then
	OUTFILTERMISMATCHNOVERREADLMAX=0.01
fi

if [ -z "$ALIGNSJSTITCHMISMATCHNMAX" ]
then
	ALIGNSJSTITCHMISMATCHNMAX="-1 -1 -1 -1"
fi


# from "params_fusions.txt":
PERC=$( cat /data/params_fusions.txt | sed -n '1 p' )				# default: 17
MINOVERLAP=$( cat /data/params_fusions.txt | sed -n '2 p' )				# default: 5
MINREADMAPPING=$( cat /data/params_fusions.txt | sed -n '3 p' )			# default: 50
STRAND=$( cat /data/params_fusions.txt | sed -n '4 p' )				# default: + (values: "+", "-" or "UN", for stranded, reverse stranded or unstranded RNASeq, respectively)

#default values
if [ -z $PERC ]
then
	PERC=17
fi

if [ -z $MINOVERLAP ]
then
	MINOVERLAP=5
fi

if [ -z $MINREADMAPPING ]
then
	MINREADMAPPING=50
fi

if [ -z $STRAND ]
then
	STRAND="+"
fi



# 3a) Collect all fusions together

echo "gene_5" "gene_3" "direct_breakpoint_5" "direct_breakpoint_3" "exon_brp_5" "exon_brp_3" "num_reads" | sed -e 's/ /\t/g' > direct_fusion_transcripts.txt

echo "gene_5" "gene_3" "direct_breakpoint_5" "direct_breakpoint_3" "exon_brp_5" "exon_brp_3" "alternative_direct_breakpoint_5" "alternative_direct_breakpoint_3" "alternative_exon_brp_5" "alternative_exon_brp_3" "num_reads" | sed -e 's/ /\t/g' > alternative_direct_fusion_transcripts.txt

echo "gene_brp_5" "gene_brp_3" "reciprocal_breakpoint_5" "reciprocal_breakpoint_3" "exon_brp_5" "exon_brp_3" "gene_bks_5" "gene_bks_3" "circular_backsplice_5" "circular_backsplice_3" "circular_exon_5" "circular_exon_3" "num_reads" | sed -e 's/ /\t/g' > reciprocal_fusion_circRNAs.txt


echo "gene_5" "gene_3" "reciprocal_breakpoint_5" "reciprocal_breakpoint_3" "exon_brp_5" "exon_brp_3" "num_reads" | sed -e 's/ /\t/g' > reciprocal_fusion_transcripts.txt

echo "gene_5" "gene_3" "reciprocal_breakpoint_5" "reciprocal_breakpoint_3" "exon_brp_5" "exon_brp_3" "alternative_reciprocal_breakpoint_5" "alternative_reciprocal_breakpoint_3" "alternative_exon_brp_5" "alternative_exon_brp_3" "num_reads" | sed -e 's/ /\t/g' > alternative_reciprocal_fusion_transcripts.txt

echo "gene_brp_5" "gene_brp_3" "direct_breakpoint_5" "direct_breakpoint_3" "exon_brp_5" "exon_brp_3" "gene_bks_5" "gene_bks_3" "circular_backsplice_5" "circular_backsplice_3" "circular_exon_5" "circular_exon_3" "num_reads" | sed -e 's/ /\t/g' > direct_fusion_circRNAs.txt



# 1) Align reads on custom reference

# run alignment
if [[ $END == 1 ]]
then

	#SE
	FILE1=$(grep $SAMPLE $SAMPLE_LIST | cut -f 2 )
	
	/tools/STAR-2.7.5c/source/STAR --runThreadN $THR --genomeDir /data/reference/genome_reference --readFilesCommand $UNCOMPRESSION_COMMAND --readFilesIn $FILE1 --outFilterMultimapNmax $OUTFILTERMULTIMAPNMAX --outFilterMatchNmin $OUTFILTERMATCHNMIN --outFilterMismatchNmax $OUTFILTERMISMATCHNMAX --outFilterMismatchNoverReadLmax $OUTFILTERMISMATCHNOVERREADLMAX --outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --alignTranscriptsPerReadNmax 20000 --alignSJstitchMismatchNmax $ALIGNSJSTITCHMISMATCHNMAX 2> "$SAMPLE"_alignment.err > "$SAMPLE"_alignment.log

else

	if [[ $END == 2 ]]
	then
	
		#PE
		FILE1=$(grep $SAMPLE $SAMPLE_LIST | cut -f 2 )
		FILE2=$(grep $SAMPLE $SAMPLE_LIST | cut -f 3 )
		
		/tools/STAR-2.7.5c/source/STAR --runThreadN $THR --genomeDir /data/reference/genome_reference --readFilesCommand $UNCOMPRESSION_COMMAND --readFilesIn $FILE1 $FILE2 --outFilterMultimapNmax $OUTFILTERMULTIMAPNMAX --outFilterMatchNmin $OUTFILTERMATCHNMIN --outFilterMismatchNmax $OUTFILTERMISMATCHNMAX --outFilterMismatchNoverReadLmax $OUTFILTERMISMATCHNOVERREADLMAX --outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --alignTranscriptsPerReadNmax 20000 --alignSJstitchMismatchNmax $ALIGNSJSTITCHMISMATCHNMAX 2> "$SAMPLE"_alignment.err > "$SAMPLE"_alignment.log

	fi

fi

# remove duplicates
samtools rmdup -s Aligned.sortedByCoord.out.bam Aligned.sortedByCoord.out_rmdup.bam

# indexing BAM file
samtools index Aligned.sortedByCoord.out_rmdup.bam

# convert BAM file in BED file
bedtools bamtobed -i Aligned.sortedByCoord.out_rmdup.bam -tag nM > Aligned.sortedByCoord.out_rmdup.bed
#convert only mapped reads



GENE1_OLD=""  #?#
GENE2_OLD=""  #?#


for ((i=1; i<=$M; i++))
do

	# retrieve fusion genes
	GENE1=$( cat /data/translocation_list2.txt | sed -n ''$i' p' | cut -f1 )
	GENE2=$( cat /data/translocation_list2.txt | sed -n ''$i' p' | cut -f2 )

	STRAND1=$( cat /data/reference/ref_$GENE1-$GENE2/info_$GENE1.txt | cut -f 5 )
	STRAND2=$( cat /data/reference/ref_$GENE1-$GENE2/info_$GENE2.txt | cut -f 5 )
	

	# retrieve the exons of the breakpoint
	EX_BRP_GENE1_1DER=$( cat /data/translocation_list2.txt | sed -n ''$i' p' | cut -f3 )
	EX_BRP_GENE2_1DER=$( cat /data/translocation_list2.txt | sed -n ''$i' p' | cut -f4 )
	EX_BRP_GENE2_2DER=$( cat /data/translocation_list2.txt | sed -n ''$i' p' | cut -f5 )
	EX_BRP_GENE1_2DER=$( cat /data/translocation_list2.txt | sed -n ''$i' p' | cut -f6 )



	# FUSIONS 1 DERIVATIVE (part 1)

	mkdir fusions_"$i"_$GENE1-$GENE2
	cd fusions_"$i"_$GENE1-$GENE2


	if [[ $GENE1 != $GENE1_OLD || $GENE2 != $GENE2_OLD ]]  #?#
	then

		# 2a) Search for fusions

		# intersection with the breakpoint junction for each gene pair
		cat ../Aligned.sortedByCoord.out_rmdup.bed | grep -e ^$GENE1"-"$GENE2"|" | sed -e 's/^.*|//' > Aligned_$GENE1-$GENE2.bed
		
		bedtools intersect -a /data/reference/ref_$GENE1-$GENE2/reference_$GENE1-$GENE2.bed -b Aligned_$GENE1-$GENE2.bed -wb > putative_fusions.bed
		
		# filter the putative fusions based on strand, the length of read mapping, the number of mismatches and the read overlap with the breakpoint
		if [[ ($STRAND == "+") || ($STRAND == "-") ]]		# stranded or reverse-stranded RNASeq
		then
			cat putative_fusions.bed | awk '$9=='\"$STRAND\"' {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $6-$5, $2-$5+1, $6-$3+1}' | sed -e 's/ /\t/g' > fusions_strand.bed
		else
			if [[ $STRAND == "UN" ]]			# unstranded RNASeq
			then
				cat putative_fusions.bed | awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $6-$5, $2-$5+1, $6-$3+1}' | sed -e 's/ /\t/g' > fusions_strand.bed
			fi
		fi

		cat fusions_strand.bed | awk '$10>='"$MINREADMAPPING"' {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $11<$12?$11:$12}' | sed -e 's/ /\t/g' > fusions_length.bed
		cat fusions_length.bed | awk '$13>='"$MINOVERLAP"' && $8<='"$PERC"'/100*$13 {print}' | sed -e 's/ /\t/g'  > fusions_filtered.bed

	else  #?#
	
		# 2a) Search for fusions
		
		cp ../fusions_"$((i-1))"_$GENE1-$GENE2/Aligned_$GENE1-$GENE2.bed .
		cp ../fusions_"$((i-1))"_$GENE1-$GENE2/putative_fusions.bed .
		cp ../fusions_"$((i-1))"_$GENE1-$GENE2/fusions_strand.bed .
		cp ../fusions_"$((i-1))"_$GENE1-$GENE2/fusions_length.bed .
		cp ../fusions_"$((i-1))"_$GENE1-$GENE2/fusions_filtered.bed .
		
	fi

	
	# 2b) Search for fusions

	# define fusion reads based on the breakpoint exons
	echo "direct_breakpoint" "mapping_breakpoint" "map_brp_junction_start" "map_brp_junction_end" "exon_map_brp_5p" "exon_map_brp_3p" "read_start" "read_end" "read_name" "n_mismatch" "library_strand" "mapping_length" "nt_left_overlap" "nt_right_overlap" | sed -e 's/ /\t/g' > direct_fusion_transcript_read.bed
	cat fusions_filtered.bed | sed -e 's/-exon/\texon/2; s/exon_//3; s/exon_//3; s/([+-])//3; s/([+-])//3' | awk '$4=='$EX_BRP_GENE1_1DER' && $5=='$EX_BRP_GENE2_1DER' {print $1, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' | sed -e 's/ /\t/g' >> direct_fusion_transcript_read.bed
	
	echo "direct_breakpoint" "mapping_breakpoint" "map_brp_junction_start" "map_brp_junction_end" "exon_map_brp_5p" "exon_map_brp_3p" "read_start" "read_end" "read_name" "n_mismatch" "library_strand" "mapping_length" "nt_left_overlap" "nt_right_overlap" > alternative_direct_fusion_transcript_read.bed
	cat fusions_filtered.bed | sed -e 's/-exon/\texon/2; s/exon_//3; s/exon_//3; s/([+-])//3; s/([+-])//3' | awk '($4<='$EX_BRP_GENE1_1DER' && $5>'$EX_BRP_GENE2_1DER') || ($4<'$EX_BRP_GENE1_1DER' && $5>='$EX_BRP_GENE2_1DER') {print "exon_'$EX_BRP_GENE1_1DER'('$STRAND1')-exon_'$EX_BRP_GENE2_1DER'('$STRAND2')", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' | sed -e 's/ /\t/g' >> alternative_direct_fusion_transcript_read.bed
	
	echo "reciprocal_breakpoint" "mapping_breakpoint" "map_brp_junction_start" "map_brp_junction_end" "exon_map_brp_5p" "exon_map_brp_3p" "read_start" "read_end" "read_name" "n_mismatch" "library_strand" "mapping_length" "nt_left_overlap" "nt_right_overlap" > reciprocal_fusion_circRNA_read.bed
	cat fusions_filtered.bed | sed -e 's/-exon/\texon/2; s/exon_//3; s/exon_//3; s/([+-])//3; s/([+-])//3' | awk '$4>='$EX_BRP_GENE1_2DER' && $5<='$EX_BRP_GENE2_2DER' {print "exon_'$EX_BRP_GENE2_2DER'('$STRAND2')-exon_'$EX_BRP_GENE1_2DER'('$STRAND1')", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' | sed -e 's/ /\t/g' >> reciprocal_fusion_circRNA_read.bed

	# summary of the number of detected reads
	echo "Number_of_fusion_transcript_reads" $( cat direct_fusion_transcript_read.bed | grep -v "breakpoint" | wc -l ) | sed -e 's/ /\t/g' > summary_reads.txt
	echo "Number_of_alternative_fusion_transcript_reads" $( cat alternative_direct_fusion_transcript_read.bed | grep -v "breakpoint" | wc -l ) | sed -e 's/ /\t/g' >> summary_reads.txt
	echo "Number_of_fusion_circRNA_reads" $( cat reciprocal_fusion_circRNA_read.bed | grep -v "breakpoint" | wc -l ) | sed -e 's/ /\t/g' >> summary_reads.txt


	# retrieve direct fusion transcripts
	cat direct_fusion_transcript_read.bed | grep -v "breakpoint" | cut -f 1,5,6 | uniq -c > file1.txt

	if [[ $( cat file1.txt | wc -l ) > 0 ]]
	then
		EX1d=$( cat file1.txt | cut -f2 )
		EX2d=$( cat file1.txt | cut -f3 )
		NRd=$( cat file1.txt | cut -f1 | sed -e 's/^ *//g; s/ /\t/' | cut -f1 )
	else
		EX1d=$EX_BRP_GENE1_1DER
		EX2d=$EX_BRP_GENE2_1DER	
	fi
	
	if [[ $STRAND1 == "+" ]]
	then
		BRP5d=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE1.bed | grep -w exon_$EX1d | awk '{print "chr"$1":"$3"("$6")"}' )
	else
		BRP5d=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE1.bed | grep -w exon_$EX1d | awk '{print "chr"$1":"$2+1"("$6")"}' )
	fi

	if [[ $STRAND2 == "+" ]]
	then
		BRP3d=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE2.bed | grep -w exon_$EX2d | awk '{print "chr"$1":"$2+1"("$6")"}' )
	else
		BRP3d=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE2.bed | grep -w exon_$EX2d | awk '{print "chr"$1":"$3"("$6")"}' )
	fi
	
	echo "gene_5" "gene_3" "direct_breakpoint_5" "direct_breakpoint_3" "exon_brp_5" "exon_brp_3" "num_reads" | sed -e 's/ /\t/g' > direct_fusion_transcript.txt

	if [[ $( cat file1.txt | wc -l ) > 0 ]]
	then
		echo $GENE1 $GENE2 $BRP5d $BRP3d $EX1d $EX2d $NRd | sed -e 's/ /\t/g' >> direct_fusion_transcript.txt
	fi


	# retrieve alternative direct fusion transcripts
	cat alternative_direct_fusion_transcript_read.bed | grep -v "breakpoint" | cut -f 1,2,5,6 | uniq -c > file2.txt
	A=$( cat file2.txt | wc -l )
	
	echo "gene_5" "gene_3" "direct_breakpoint_5" "direct_breakpoint_3" "exon_brp_5" "exon_brp_3" "alternative_direct_breakpoint_5" "alternative_direct_breakpoint_3" "alternative_exon_brp_5" "alternative_exon_brp_3" "num_reads" | sed -e 's/ /\t/g' > alternative_direct_fusion_transcript.txt

	for ((a=1; a<=$A; a++))
	do
		EX1a=$( cat file2.txt | sed -n $a'p' | cut -f3 )
		EX2a=$( cat file2.txt | sed -n $a'p' | cut -f4 )
		NRa=$( cat file2.txt | sed -n $a'p' | cut -f1 | sed -e 's/^ *//g; s/ /\t/' | cut -f1 )

		if [[ $STRAND1 == "+" ]]
		then
			BRP5a=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE1.bed | grep -w exon_$EX1a | awk '{print "chr"$1":"$3"("$6")"}' )
		else
			BRP5a=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE1.bed | grep -w exon_$EX1a | awk '{print "chr"$1":"$2+1"("$6")"}' )
		fi

		if [[ $STRAND2 == "+" ]]
		then
			BRP3a=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE2.bed | grep -w exon_$EX2a | awk '{print "chr"$1":"$2+1"("$6")"}' )
		else
			BRP3a=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE2.bed | grep -w exon_$EX2a | awk '{print "chr"$1":"$3"("$6")"}' )
		fi
		
		echo $GENE1 $GENE2 $BRP5d $BRP3d $EX1d $EX2d $BRP5a $BRP3a $EX1a $EX2a $NRa | sed -e 's/ /\t/g' >> alternative_direct_fusion_transcript.txt
	done



	cd ..



	# FUSIONS 2 DERIVATIVE

	mkdir fusions_"$i"_$GENE2-$GENE1
	cd fusions_"$i"_$GENE2-$GENE1


	if [[ $GENE1 != $GENE1_OLD || $GENE2 != $GENE2_OLD ]]  #?#
	then

		# 2a) Search for fusions

		# intersection with the backsplice junction for each gene pair
		cat ../Aligned.sortedByCoord.out_rmdup.bed | grep -e ^$GENE2"-"$GENE1"|" | sed -e 's/^.*|//' > Aligned_$GENE2-$GENE1.bed
		
		bedtools intersect -a /data/reference/ref_$GENE1-$GENE2/reference_$GENE2-$GENE1.bed -b Aligned_$GENE2-$GENE1.bed -wb > putative_fusions.bed

		# filter the putative fusions based on strand, the length of read mapping, the number of mismatch and the read overlap with the breakpoint
		if [[ ($STRAND == "+") || ($STRAND == "-") ]]		# stranded or reverse-stranded RNASeq
		then
			cat putative_fusions.bed | awk '$9=='\"$STRAND\"' {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $6-$5, $2-$5+1, $6-$3+1}' | sed -e 's/ /\t/g' > fusions_strand.bed
		else
			if [[ $STRAND == "UN" ]]			# unstranded RNASeq
			then
				cat putative_fusions.bed | awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $6-$5, $2-$5+1, $6-$3+1}' | sed -e 's/ /\t/g' > fusions_strand.bed
			fi
		fi

		cat fusions_strand.bed | awk '$10>='"$MINREADMAPPING"' {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $11<$12?$11:$12}' | sed -e 's/ /\t/g' > fusions_length.bed
		cat fusions_length.bed | awk '$13>='"$MINOVERLAP"' && $8<='"$PERC"'/100*$13 {print}' | sed -e 's/ /\t/g'  > fusions_filtered.bed

	else  #?#

		# 2a) Search for fusions

		cp ../fusions_"$((i-1))"_$GENE2-$GENE1/Aligned_$GENE2-$GENE1.bed .
		cp ../fusions_"$((i-1))"_$GENE2-$GENE1/putative_fusions.bed .
		cp ../fusions_"$((i-1))"_$GENE2-$GENE1/fusions_strand.bed .
		cp ../fusions_"$((i-1))"_$GENE2-$GENE1/fusions_length.bed .
		cp ../fusions_"$((i-1))"_$GENE2-$GENE1/fusions_filtered.bed .
		
	fi
	
			
	# 2b) Search for fusions
	
	# define fusion reads based on the breakpoint exons
	echo "reciprocal_breakpoint" "mapping_breakpoint" "map_brp_junction_start" "map_brp_junction_end" "exon_map_brp_5p" "exon_map_brp_3p" "read_start" "read_end" "read_name" "n_mismatch" "library_strand" "mapping_length" "nt_left_overlap" "nt_right_overlap" | sed -e 's/ /\t/g' > reciprocal_fusion_transcript_read.bed
	cat fusions_filtered.bed | sed -e 's/-exon/\texon/2; s/exon_//3; s/exon_//3; s/([+-])//3; s/([+-])//3' | awk '$4=='$EX_BRP_GENE2_2DER' && $5=='$EX_BRP_GENE1_2DER' {print $1, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' | sed -e 's/ /\t/g' >> reciprocal_fusion_transcript_read.bed

	echo "reciprocal_breakpoint" "mapping_breakpoint" "map_brp_junction_start" "map_brp_junction_end" "exon_map_brp_5p" "exon_map_brp_3p" "read_start" "read_end" "read_name" "n_mismatch" "library_strand" "mapping_length" "nt_left_overlap" "nt_right_overlap" > alternative_reciprocal_fusion_transcript_read.bed
	cat fusions_filtered.bed | sed -e 's/-exon/\texon/2; s/exon_//3; s/exon_//3; s/([+-])//3; s/([+-])//3' | awk '($4<='$EX_BRP_GENE2_2DER' && $5>'$EX_BRP_GENE1_2DER') || ($4<'$EX_BRP_GENE2_2DER' && $5>='$EX_BRP_GENE1_2DER') {print "exon_'$EX_BRP_GENE2_2DER'('$STRAND2')-exon_'$EX_BRP_GENE1_2DER'('$STRAND1')", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' | sed -e 's/ /\t/g' >> alternative_reciprocal_fusion_transcript_read.bed

	echo "direct_breakpoint" "mapping_breakpoint" "map_brp_junction_start" "map_brp_junction_end" "exon_map_brp_5p" "exon_map_brp_3p" "read_start" "read_end" "read_name" "n_mismatch" "library_strand" "mapping_length" "nt_left_overlap" "nt_right_overlap" > direct_fusion_circRNA_read.bed
	cat fusions_filtered.bed | sed -e 's/-exon/\texon/2; s/exon_//3; s/exon_//3; s/([+-])//3; s/([+-])//3' | awk '$4>='$EX_BRP_GENE2_1DER' && $5<='$EX_BRP_GENE1_1DER' {print "exon_'$EX_BRP_GENE1_1DER'('$STRAND1')-exon_'$EX_BRP_GENE2_1DER'('$STRAND2')", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' | sed -e 's/ /\t/g' >> direct_fusion_circRNA_read.bed

	# summary of the number of detected reads
	echo "Number_of_fusion_transcript_reads" $( cat reciprocal_fusion_transcript_read.bed | grep -v "breakpoint" | wc -l ) | sed -e 's/ /\t/g' > summary_reads.txt
	echo "Number_of_alternative_fusion_transcript_reads" $( cat alternative_reciprocal_fusion_transcript_read.bed | grep -v "breakpoint" | wc -l ) | sed -e 's/ /\t/g' >> summary_reads.txt
	echo "Number_of_fusion_circRNA_reads" $( cat direct_fusion_circRNA_read.bed | grep -v "breakpoint" | wc -l ) | sed -e 's/ /\t/g' >> summary_reads.txt


	# retrieve reciprocal fusion transcripts
	cat reciprocal_fusion_transcript_read.bed| grep -v "breakpoint" | cut -f 1,5,6 | uniq -c > file1.txt

	if [[ $( cat file1.txt | wc -l ) > 0 ]]
	then
		EX1r=$( cat file1.txt | cut -f2 )
		EX2r=$( cat file1.txt | cut -f3 )
		NRr=$( cat file1.txt | cut -f1 | sed -e 's/^ *//g; s/ /\t/' | cut -f1 )
	else
		EX1r=$EX_BRP_GENE2_2DER
		EX2r=$EX_BRP_GENE1_2DER
	fi	
	
	if [[ $STRAND2 == "+" ]]
	then
		BRP5r=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE2.bed | grep -w exon_$EX1r | awk '{print "chr"$1":"$3"("$6")"}' )
	else
		BRP5r=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE2.bed | grep -w exon_$EX1r | awk '{print "chr"$1":"$2+1"("$6")"}' )
	fi

	if [[ $STRAND1 == "+" ]]
	then
		BRP3r=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE1.bed | grep -w exon_$EX2r | awk '{print "chr"$1":"$2+1"("$6")"}' )
	else
		BRP3r=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE1.bed | grep -w exon_$EX2r | awk '{print "chr"$1":"$3"("$6")"}' )
	fi
	
	echo "gene_5" "gene_3" "reciprocal_breakpoint_5" "reciprocal_breakpoint_3" "exon_brp_5" "exon_brp_3" "num_reads" | sed -e 's/ /\t/g' > reciprocal_fusion_transcript.txt
	
	if [[ $( cat file1.txt | wc -l ) > 0 ]]
	then
		echo $GENE2 $GENE1 $BRP5r $BRP3r $EX1r $EX2r $NRr | sed -e 's/ /\t/g' >> reciprocal_fusion_transcript.txt
	fi


	# retrieve alternative reciprocal fusion transcripts
	cat alternative_reciprocal_fusion_transcript_read.bed | grep -v "breakpoint" | cut -f 1,2,5,6 | uniq -c > file2.txt
	A=$( cat file2.txt | wc -l )
	
	echo "gene_5" "gene_3" "reciprocal_breakpoint_5" "reciprocal_breakpoint_3" "exon_brp_5" "exon_brp_3" "alternative_reciprocal_breakpoint_5" "alternative_reciprocal_breakpoint_3" "alternative_exon_brp_5" "alternative_exon_brp_3" "num_reads" | sed -e 's/ /\t/g' > alternative_reciprocal_fusion_transcript.txt

	for ((a=1; a<=$A; a++))
	do
		EX1a=$( cat file2.txt | sed -n $a'p' | cut -f3 )
		EX2a=$( cat file2.txt | sed -n $a'p' | cut -f4 )
		NRa=$( cat file2.txt | sed -n $a'p' | cut -f1 | sed -e 's/^ *//g; s/ /\t/' | cut -f1 )

		if [[ $STRAND2 == "+" ]]
		then
			BRP5a=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE2.bed | grep -w exon_$EX1a | awk '{print "chr"$1":"$3"("$6")"}' )
		else
			BRP5a=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE2.bed | grep -w exon_$EX1a | awk '{print "chr"$1":"$2+1"("$6")"}' )
		fi

		if [[ $STRAND1 == "+" ]]
		then
			BRP3a=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE1.bed | grep -w exon_$EX2a | awk '{print "chr"$1":"$2+1"("$6")"}' )
		else
			BRP3a=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE1.bed | grep -w exon_$EX2a | awk '{print "chr"$1":"$3"("$6")"}' )
		fi
		
		echo $GENE2 $GENE1 $BRP5r $BRP3r $EX1r $EX2r $BRP5a $BRP3a $EX1a $EX2a $NRa | sed -e 's/ /\t/g' >> alternative_reciprocal_fusion_transcript.txt
	done


	# retrieve direct fusion circRNAs
	cat direct_fusion_circRNA_read.bed | grep -v "breakpoint" | cut -f 1,2,5,6 | uniq -c > file3.txt
	C=$( cat file3.txt | wc -l )
	
	echo "gene_brp_5" "gene_brp_3" "direct_breakpoint_5" "direct_breakpoint_3" "exon_brp_5" "exon_brp_3" "gene_bks_5" "gene_bks_3" "circular_backsplice_5" "circular_backsplice_3" "circular_exon_5" "circular_exon_3" "num_reads" | sed -e 's/ /\t/g' > direct_fusion_circRNA.txt

	for ((c=1; c<=$C; c++))
	do
		EX1c=$( cat file3.txt | sed -n $c'p' | cut -f3 )
		EX2c=$( cat file3.txt | sed -n $c'p' | cut -f4 )
		NRc=$( cat file3.txt | sed -n $c'p' | cut -f1 | sed -e 's/^ *//g; s/ /\t/' | cut -f1 )

		if [[ $STRAND2 == "+" ]]
		then
			BRP5c=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE2.bed | grep -w exon_$EX1c | awk '{print "chr"$1":"$3"("$6")"}' )
		else
			BRP5c=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE2.bed | grep -w exon_$EX1c | awk '{print "chr"$1":"$2+1"("$6")"}' )
		fi

		if [[ $STRAND1 == "+" ]]
		then
			BRP3c=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE1.bed | grep -w exon_$EX2c | awk '{print "chr"$1":"$2+1"("$6")"}' )
		else
			BRP3c=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE1.bed | grep -w exon_$EX2c | awk '{print "chr"$1":"$3"("$6")"}' )
		fi
		
		echo $GENE1 $GENE2 $BRP5d $BRP3d $EX1d $EX2d $GENE2 $GENE1 $BRP5c $BRP3c $EX1c $EX2c $NRc | sed -e 's/ /\t/g' >> direct_fusion_circRNA.txt
	done



	# 3b) Collect all fusions together

	cat reciprocal_fusion_transcript.txt | grep -v "breakpoint" >> ../reciprocal_fusion_transcripts_part.txt
	cat alternative_reciprocal_fusion_transcript.txt | grep -v "breakpoint" >> ../alternative_reciprocal_fusion_transcripts_part.txt
	cat direct_fusion_circRNA.txt | grep -v "breakpoint" >> ../direct_fusion_circRNAs_part.txt

	# summary of the number of detected fusions
	echo "Number_of_fusion_transcripts" $( cat reciprocal_fusion_transcript.txt | grep -v "breakpoint" | wc -l ) | sed -e 's/ /\t/g' > summary_fusions.txt
	echo "Number_of_alternative_fusion_transcripts" $( cat alternative_reciprocal_fusion_transcript.txt | grep -v "breakpoint" | wc -l ) | sed -e 's/ /\t/g' >> summary_fusions.txt
	echo "Number_of_fusion_circRNAs" $( cat direct_fusion_circRNA.txt | grep -v "breakpoint" | wc -l ) | sed -e 's/ /\t/g' >> summary_fusions.txt


	rm file1.txt
	rm file2.txt
	rm file3.txt



	cd ..



	# FUSIONS 1 DERIVATIVE (part 2)

	cd fusions_"$i"_$GENE1-$GENE2

	# 2b) Search for fusions

	# retrieve reciprocal fusion circRNAs
	cat reciprocal_fusion_circRNA_read.bed | grep -v "breakpoint" | cut -f 1,2,5,6 | uniq -c > file3.txt
	C=$( cat file3.txt | wc -l )
	
	echo "gene_brp_5" "gene_brp_3" "reciprocal_breakpoint_5" "reciprocal_breakpoint_3" "exon_brp_5" "exon_brp_3" "gene_bks_5" "gene_bks_3" "circular_backsplice_5" "circular_backsplice_3" "circular_exon_5" "circular_exon_3" "num_reads" | sed -e 's/ /\t/g' > reciprocal_fusion_circRNA.txt

	for ((c=1; c<=$C; c++))
	do
		EX1c=$( cat file3.txt | sed -n $c'p' | cut -f3 )
		EX2c=$( cat file3.txt | sed -n $c'p' | cut -f4 )
		NRc=$( cat file3.txt | sed -n $c'p' | cut -f1 | sed -e 's/^ *//g; s/ /\t/' | cut -f1 )

		if [[ $STRAND1 == "+" ]]
		then
			BRP5c=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE1.bed | grep -w exon_$EX1c | awk '{print "chr"$1":"$3"("$6")"}' )
		else
			BRP5c=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE1.bed | grep -w exon_$EX1c | awk '{print "chr"$1":"$2+1"("$6")"}' )
		fi

		if [[ $STRAND2 == "+" ]]
		then
			BRP3c=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE2.bed | grep -w exon_$EX2c | awk '{print "chr"$1":"$2+1"("$6")"}' )
		else
			BRP3c=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE2.bed | grep -w exon_$EX2c | awk '{print "chr"$1":"$3"("$6")"}' )
		fi
		
		echo $GENE2 $GENE1 $BRP5r $BRP3r $EX1r $EX2r $GENE1 $GENE2 $BRP5c $BRP3c $EX1c $EX2c $NRc | sed -e 's/ /\t/g' >> reciprocal_fusion_circRNA.txt
	done



	# 3b) Collect all fusions together

	cat direct_fusion_transcript.txt | grep -v "breakpoint" >> ../direct_fusion_transcripts_part.txt
	cat alternative_direct_fusion_transcript.txt | grep -v "breakpoint" >> ../alternative_direct_fusion_transcripts_part.txt
	cat reciprocal_fusion_circRNA.txt | grep -v "breakpoint" >> ../reciprocal_fusion_circRNAs_part.txt

	# summary of the number of detected fusions
	echo "Number_of_fusion_transcripts" $( cat direct_fusion_transcript.txt | grep -v "breakpoint" | wc -l ) | sed -e 's/ /\t/g' > summary_fusions.txt
	echo "Number_of_alternative_fusion_transcripts" $( cat alternative_direct_fusion_transcript.txt | grep -v "breakpoint" | wc -l ) | sed -e 's/ /\t/g' >> summary_fusions.txt
	echo "Number_of_fusion_circRNAs" $( cat reciprocal_fusion_circRNA.txt | grep -v "breakpoint" | wc -l ) | sed -e 's/ /\t/g' >> summary_fusions.txt


	rm file1.txt
	rm file2.txt
	rm file3.txt
	
	

	cd ..
	
	

GENE1_OLD=$GENE1  #?#
GENE2_OLD=$GENE2



done



# 3c) Collect all fusions together

# remove duplicated fusions
cat direct_fusion_transcripts_part.txt | sort | uniq >> direct_fusion_transcripts.txt
cat alternative_direct_fusion_transcripts_part.txt | sort | uniq >> alternative_direct_fusion_transcripts.txt
cat reciprocal_fusion_circRNAs_part.txt | sort | uniq >> reciprocal_fusion_circRNAs.txt

cat reciprocal_fusion_transcripts_part.txt | sort | uniq >> reciprocal_fusion_transcripts.txt
cat alternative_reciprocal_fusion_transcripts_part.txt | sort | uniq >> alternative_reciprocal_fusion_transcripts.txt
cat direct_fusion_circRNAs_part.txt | sort | uniq >> direct_fusion_circRNAs.txt


# remove unuseless files
rm direct_fusion_transcripts_part.txt
rm alternative_direct_fusion_transcripts_part.txt
rm reciprocal_fusion_circRNAs_part.txt

rm reciprocal_fusion_transcripts_part.txt
rm alternative_reciprocal_fusion_transcripts_part.txt
rm direct_fusion_circRNAs_part.txt


# create .csv files
LIST=$(ls *.txt)

for FILE in $LIST
do
	cat $FILE | sed -e 's/\t/,/g' > ${FILE%.txt}.csv
done

