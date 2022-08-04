# CircFusion

CircFusion is a new software tool to detect linear fusions and f-circRNAs from RNA-seq data, both in samples where the fusion breakpoint is known and in a discovery mode. CircFusion can detect chimeric transcripts even in the presence of multiple translocation breakpoints, both for the main and unbalanced reciprocal fusion, which are quite common in malignant cells. The CircFusion output includes an HTML report for each sample, with the f-circRNAs and linear RNAs identified, summarized for each translocation. Plus, it gives plain text files with the genomic coordinates of the identified fusions, and information about the read alignment. 

## Installation

### Installation from the Docker image

The Docker image saves you from the installation burden. A Docker image of CircFusion is available from DockerHub; just pull it with the command:

    docker pull annadalmolin/circfusion:v1.0
    
## Usage

### Input data

Prepare your project directory with the following files:

- _translocation_list.txt_: file with the list of translocations to be searched. The file format is a tab-separated text file, with the following 8 fields: 
  1. name of the gene in 5’ in the translocation (the official gene name has to be used)
  2. name of the gene in 3’ in the translocation (the official gene name has to be used)
  3. exon number of the gene in 5’ where the direct breakpoint occurs (optional)
  4. exon number of the gene in 3’ where the direct breakpoint occurs (optional)
  5. exon number of the gene in 3’ where the reciprocal breakpoint occurs (optional)
  6. exon number of the gene in 5’ where the reciprocal breakpoint occurs (optional)
  7. Ensembl id according to the transcript of the gene in 5’
  8. Ensembl id according to the transcript of the gene in 3’

  Tips:
  - when the exon number is not known, write a dot in the corresponding field of the file
  - gene pairs fusing the same two genes in multiple breakpoints should be listed one after the other

  
  An example of _translocation_list.txt_ is:
  
      BCR      ABL1      14   2    1    15    ENST00000305877    ENST00000372348
      KMT2A    MLLT1     9    2    .    .     ENST00000534358    ENST00000252674
      KMT2A    MLLT10    .    .    .    .     ENST00000534358    ENST00000377072
      SLC34A2  ROS1      4    34   33   5     ENST00000382051    ENST00000368507
      SLC34A2  ROS1      4    32   31   5     ENST00000382051    ENST00000368507
      
  Note: the reciprocal fusion exons are indicated in the same line of the direct fusion exons for the same gene pair, and not in a new line.

- _sample_list.txt_: file with the list of samples to be analyzed. The file format is a tab-separated text file, with the sample name in the first column and the path to the RNA-seq file in the second (for single-end data) or in the second and third columns (for paired-end data). Sample order doesn’t matter.

  An example of _sample_list.txt_ (SE) is:

      SAMPLE1    /data/reads/sample1_SE.fq.gz
      SAMPLE2    /data/reads/sample2_SE.fq.gz

  An example of _sample_list.txt_ (PE) is:

      SAMPLE1           /data/reads/sample1_PE_1.fq.gz        /data/reads/sample1_PE_2.fq.gz
      SAMPLE2           /data/reads/sample2_PE_1.fq.gz        /data/reads/sample2_PE_2.fq.gz
     
- _path_file.txt_: file with the paths for Ensembl annotation and genome files. The file format is a text file with a path written in each row, __in the following order__:
  1. path to annotation file
  2. path to genome file
  
  An example of _path_files.txt_ is:
  
      /data/annotation/Homo_sapiens.GRCh38.104.gtf
      /data/annotation/Homo_sapiens.GRCh38.dna.primary_assembly.fa

  The gene annotation (in GTF format) and the genome sequence (in FASTA format) files must be downloaded by the user from the Ensembl database and placed into the _annotation_/ directory contained in the project directory. Annotation and genome files for Homo sapiens (GRCh38) can be downloaded from http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/ and http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/, respectively.

- _params_Star.txt_: file with the parameters used by Star embedded in CircFusion. The file format is a text file with a parameter written in each row, __in the following order__:

  1. number of threads for the analysis;
  2. uncompression command; the default is _“zcat”_ for gzip compressed file (.gz);
  3. maximum number of multiple alignments allowed (default value 10);
  4. minimum number of matches required for the alignment (default value 5);
  5. maximum number of mismatches allowed in the alignment (default value 4);
  6. maximum number of mismatches allowed per pair relative to read length (default value 0.01);
  7. maximum number of mismatches allowed for stitching of the splice junctions (default value -1 -1 -1 -1).

  An example of _params_Star.txt_ is:
  
      10
      "zcat"
      10
       
       
      0.01
      -1 -1 -1 -1

- _params_fusions.txt_: file with the parameters used by CircFusion to define fusion transcripts. The file format is a text file with a parameter written in each row, __in the following order__:

  1. maximum percentage of mismatches allowed at the size of shorter overlap over the fusion junction (default value 17);
  2. minimum number of nucleotides overlapping the fusion junction (default value 5);
  3. minimum number of nucleotides mapping in the alignment (default value 50);
  4. library strand; allowed values are "+", "-" or "UN", for stranded, reverse stranded or unstranded RNASeq, respectively (default value “+”).

  An example of _params_fusions.txt_ is:

      17
      5
       
      +

and directories:

- _annotation/_: directory containing the following files:

  1. gene annotation file from the Ensembl database
  2. genome sequence file from the Ensembl database
  
- _reads/_: directory containing the RNA-seq files to be aligned. Both FASTA and FASTQ files are allowed. If the read files are compressed, specify the uncompression command in the _params_Star.txt_ file; different compression formats are allowed, according to Star aligner  (https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).


### Missing information about the breakpoint exons

If the _translocation_list.txt_ file has full information about the breakpoint exons, i.e. the user knows the exons where the direct and reciprocal breakpoints occur for a gene pair, CircFusion will use this information to define fusion transcripts. 
If some information is missing, CircFusion algorithm retrieves it automatically:

- when only the exons of the direct translocation are indicated, the exons of the reciprocal translocation are retrieved considering a corresponding reciprocal translocation compatible with the fusion;
- when only the gene pair is specified in the file, without information about the fused exons, then the complete reference sequences are scanned, and CircFusion treats in turn each pair of exons of those genes as the putative direct and reciprocal exons, defining the fusion transcripts compatible with those breakpoints.


### Running the analysis

To run CircFusion from the Docker container use:

    sudo docker run -it -v $(pwd):/data annadalmolin/circfusion:v1.0
  
All paths in _sample_list.txt_ and _path_files.txt_ must be relative to the directory in the container where the volumes were mounted (f.i. /data/reads/sample1_SE.fq.gz, as detailed above). 

If you want the container to give your user permissions, you need to set the owner id with "-u `id -u`":

    sudo docker run -u `id -u` -it -v $(pwd):/data annadalmolin/circfusion:v1.0


### Output data

After CircFusion successful run end, you will find the following new directories in your project directory:
1. _reference/_: contains the reference sequences created by CircFusion and used for the alignment
2. _fusions/_: contains the detected linear fusions and f-circRNAs
3. _graphical_output/_: contains the HTML summary files with the detected fusions

- __reference/__

  The directory contains:
  - _reference_sequences.fa_: FASTA file with the reference sequences generated by exon combination for all gene pairs in _translocation_list.txt_
  - _genome_reference/_: directory with the Star index
  - as many subdirectories as the unique pairs of fused genes in the _translocation_list.txt_ file. If two genes are fused by more breakpoints (both direct and reciprocals), the reference sequences are included in the same subdirectory.
Each subdirectory contains:
    - text and BED files about the genomic coordinates and the exon composition of the 5’ gene of the translocation
    - text and BED files about the genomic coordinates and the exon composition of the 3’ gene of the translocation
    - text, BED and FASTA files about the exon relative position and exon sequences of both the 5’ gene fused with the 3’ gene, and the 3’ gene fused with the 5’ gene of the translocation

- __fusions/__

  The directory contains as many subdirectories as the samples investigated. Each sample subdirectory includes:
  - separated text files for the detected fusion RNAs, considering all translocations together (start positions are 1-based):
    - _direct_fusion_transcripts.txt_ and _reciprocal_fusion_transcripts.txt_: linear fusions from the direct or the reciprocal translocation, respectively;
    - _alternative_direct_fusion_transcripts.txt_ and _alternative_reciprocal_fusion_transcripts.txt_: linear fusions from the direct or the reciprocal translocation, respectively, according to breakpoints different from those indicated by the user;
    - _reciprocal_fusion_circRNAs.txt_ and _direct_fusion_circRNAs.txt_: fusion circRNAs according to the direct or the reciprocal translocation, respectively;

  - two subdirectories for each investigated translocation, one for the direct and one for the reciprocal translocation, with separated text files for the detected fusion RNAs and reads, according to that translocation:
    - files with the putative fusion reads filtered step by step according to the parameters in _params_fusions.txt_ (start positions are 0-based):
      - _putative_fusions.bed_: all fusion reads aligned to the reference and crossing the fusion junction;
      - _fusions_strand.bed_: reads from _putative_fusions.bed_ filtered for the library strand;
      - _fusions_length.bed_: reads from _fusions_strand.bed_ with enough nucleotides mapping in the alignment;
      - _fusions_filtered.bed_: reads from _fusions_length.bed_ filtered for percentage of mismatches and number of nucleotides overlapping the fusion junction;

    - files with the selected fusion reads according to the breakpoint exons (start positions are 0-based):
      - _direct_fusion_transcript_read.bed_ and _reciprocal_fusion_transcript_read.bed_: supporting reads for the linear fusions from the direct or the reciprocal translocation, respectively;
      - _alternative_direct_fusion_transcript_read.bed_ and _alternative_reciprocal_fusion_transcript_read.bed_: supporting reads for the linear fusions from the direct or the reciprocal translocation, respectively, according to breakpoints different from those indicated by the user;
      - _reciprocal_fusion_circRNA_read.bed_ and _direct_fusion_circRNA_read.bed_: supporting reads for the fusion circRNAs according to the direct or the reciprocal translocation, respectively;

    - files with the fusion RNAs defined by the fusion reads:
      - _direct_fusion_transcript.txt_ and _reciprocal_fusion_transcript.txt_: linear fusions from the direct or the reciprocal translocation, respectively;
      - _alternative_direct_fusion_transcript.txt_ and _alternative_reciprocal_fusion_transcript.txt_: linear fusions from the direct or the reciprocal translocation, respectively, according to breakpoints different from those indicated by the user;
      - _reciprocal_fusion_circRNA.txt_ and _direct_fusion_circRNA.txt_: fusion circRNAs according to the direct or the reciprocal translocation, respectively;

    - two summary text files:
      - _summary_reads.txt_: counting the different fusion reads;
      - _summary_fusions.txt_: counting the different fusion transcripts;

  - output files of Star aligner:
      - log and error files;
      - files with aligned reads in BAM and BED formats;
      - file with unmapped reads.

- __graphical_output/__

  The directory contains as many subdirectories as the samples investigated. Each sample subdirectory includes:
  - an HTML file with all the detected fusions, separated for type and translocation
  - separated plots for each type of fusion detected (fusion transcripts, alternative fusion transcripts, or f-circRNAs; direct or reciprocal)


## Additional notes

- The “alternative fusion transcripts” category denotes additional translocations in the sample. In this case, the user can run the analysis again, indicating these newly discovered translocations in the _translocation_list.txt_ file, to search for further f-circRNAs and linear fusions compatible with the new junctions.
- Before running the analysis again on the same dataset:
  - remove the _translocation_list2.txt_ file
  - keep the _reference/_ directory if you want to search for the same fused gene pair
  - remove or rename the _fusion/_ directory
  - remove or rename the _graphical_output/_ directory
