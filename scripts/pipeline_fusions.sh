#!/bin/bash

## Main script for linear and circular fusion detection

## input files:
# file with the list of translocations to be searched: "translocation_list.txt"
# file with the list of samples to be analyzed: "sample_list.txt"
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
# 
#"sample_list.txt":
#1 field: sample name
#2 field: path to SE data or to PE_1 data
#3 field: path to PE_2 data (optional)

## command: ./scripts/pipeline_fusions.sh translocation_list.txt sample_list.txt



TRANSLOCATIONS=$1
SAMPLE_LIST=$2

THR=$( cat /data/params_Star.txt | sed -n '1 p' )					# default: 8

#default value
if [ -z $THR ]
	then
		THR=8
	fi



# 1) Create the reference file for each translocation

if [ ! -d "/data/reference" ]
then
	mkdir /data/reference
	cd /data/reference

	# generate reference sequences	
	cat $TRANSLOCATIONS | cut -f1,2 | uniq > trans.txt

	N=$( cat trans.txt | wc -l )

	for ((i=1; i<=$N; i++))
	do

		GENE1=$( cat trans.txt | sed -n ''$i' p' | cut -f1 )
		GENE2=$( cat trans.txt | sed -n ''$i' p' | cut -f2 )

		mkdir ref_$GENE1-$GENE2
		cd ref_$GENE1-$GENE2

		/scripts/script_find_fusions_Star_part1.sh $GENE1 $GENE2
		
		cd ..

	done


	# generate STAR index (small genomes):
	/tools/STAR-2.7.5c/source/STAR --runMode genomeGenerate --runThreadN $THR --genomeDir genome_reference --genomeFastaFiles reference_sequences.fa --genomeSAindexNbases 10


	rm trans.txt

	cd ..

else
	echo "Using the reference sequence set already present."
fi



# 2) Detect fusions for each sample

if [ ! -d "/data/fusions" ]
then

	mkdir /data/fusions
	cd /data/fusions


	# number of translocation to detect
	M=$( cat $TRANSLOCATIONS | wc -l )

	for ((i=1; i<=$M; i++))
	do

		# retrieve fusion genes
		GENE1=$( cat $TRANSLOCATIONS | sed -n ''$i' p' | cut -f1 )
		GENE2=$( cat $TRANSLOCATIONS | sed -n ''$i' p' | cut -f2 )	

		ENS1=$( cat $TRANSLOCATIONS | sed -n ''$i' p' | cut -f7 )
		ENS2=$( cat $TRANSLOCATIONS | sed -n ''$i' p' | cut -f8 )

		# retrieve the exons of the breakpoint
		EX_BRP_GENE1_1DER=$( cat $TRANSLOCATIONS | sed -n ''$i' p' | cut -f3 )
		EX_BRP_GENE2_1DER=$( cat $TRANSLOCATIONS | sed -n ''$i' p' | cut -f4 )
		EX_BRP_GENE2_2DER=$( cat $TRANSLOCATIONS | sed -n ''$i' p' | cut -f5 )
		EX_BRP_GENE1_2DER=$( cat $TRANSLOCATIONS | sed -n ''$i' p' | cut -f6 )

		if [[ ($EX_BRP_GENE1_1DER != ".") && ($EX_BRP_GENE2_1DER != ".") ]]	# if exist
		then
			if [ $EX_BRP_GENE2_2DER == "." ]				# if doesn't exist
			then
				EX_BRP_GENE2_2DER=$(($EX_BRP_GENE2_1DER-1))		# balanced reciprocal translocation
			fi

			if [ $EX_BRP_GENE1_2DER == "." ]
			then
				EX_BRP_GENE1_2DER=$(($EX_BRP_GENE1_1DER+1))		# balanced reciprocal translocation
			fi
			
			echo $GENE1 $GENE2 $EX_BRP_GENE1_1DER $EX_BRP_GENE2_1DER $EX_BRP_GENE2_2DER $EX_BRP_GENE1_2DER $ENS1 $ENS2 | sed -e 's/ /\t/g' >> ../translocation_list2.txt
			
		else

			# make the combinations for each exon pair for the two genes
			N_EX_GENE1=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE1.bed | wc -l )
			N_EX_GENE2=$( cat /data/reference/ref_$GENE1-$GENE2/exons_$GENE2.bed | wc -l )
		
			for ((G1=1; G1<=$(($N_EX_GENE1-1)); G1++))
			do
				EX_BRP_GENE1_1DER=$G1
				
				for ((G2=2; G2<=$N_EX_GENE2; G2++))
				do
					EX_BRP_GENE2_1DER=$G2
					EX_BRP_GENE2_2DER=$(($EX_BRP_GENE2_1DER-1))		# balanced reciprocal translocation
					EX_BRP_GENE1_2DER=$(($EX_BRP_GENE1_1DER+1))		# balanced reciprocal translocation

					echo $GENE1 $GENE2 $EX_BRP_GENE1_1DER $EX_BRP_GENE2_1DER $EX_BRP_GENE2_2DER $EX_BRP_GENE1_2DER $ENS1 $ENS2 | sed -e 's/ /\t/g' >> ../translocation_list2.txt

				done
				
			done
							
		fi

	done

	M=$( cat /data/translocation_list2.txt | wc -l )


	# SE or PE reads
	NFS=$( cat $SAMPLE_LIST | awk -F "\t" '{print NF}' | uniq )
	END=$(($NFS-1))


	if [[ ($END == 1) || ($END == 2) ]]
	then

		for SAMPLE in $( cat $SAMPLE_LIST | cut -f 1 )
		do

			mkdir $SAMPLE
			cd $SAMPLE
			
			/scripts/script_find_fusions_Star_part2.sh $M $END $SAMPLE $SAMPLE_LIST
			
			cd ..
			
		done
		
	else

		echo "Uncorrect column number in sample_list.txt file"

	fi


	cd ..

else

	echo "Remove or rename the \"fusion\" directory before running the analysis again!"

fi



# 3) Graphical output

# prepare 1
if [ ! -d "/data/.R" ]
then
	mkdir /data/.R
fi

if [ ! -d "/data/.R/lib" ]
then
	mkdir /data/.R/lib
fi


# prepare 2
if [ ! -d "/data/.scripts" ]
then
	mkdir /data/.scripts
fi

cp /scripts/script_detected_fusions.Rmd /data/.scripts/



if [ ! -d "/data/graphical_output" ]
then

	mkdir /data/graphical_output
	cd /data/graphical_output

	for SAMPLE in $( cat $SAMPLE_LIST | cut -f 1 )
	do
		mkdir $SAMPLE
		cd $SAMPLE

		R -e "rmarkdown::render('/data/.scripts/script_detected_fusions.Rmd', output_file='script_detected_fusions.html', output_dir='.', params = list(sample='$SAMPLE'))"
		
		cd ..

	done

	cd ..

else

	echo "Remove or rename the \"graphical_output\" directory before running the analysis again!"

fi



# remove unuseless files
rm -rf /data/.scripts
rm -rf /data/.R

