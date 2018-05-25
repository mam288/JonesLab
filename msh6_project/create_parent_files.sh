#!/bin/bash

# Input Parameters
steps=$1 # which processing steps to perform Align: Sort: Extract: Sort and Index: Create Coverage Files

# Set variables
# CHANGE THESE IF DESIRED
ct=5 #coverage threshold
ad=2 #allele depth threshold
mq=0 #mapping quality threshold
par1_name="par1" # name of the first parent fly line
par2_name="par2" # name of the second parent fly line
#offspring_name="offspring" # name of the offspring fly line
ref="dm6" # ref genome name

# File Paths
# CHANGE THESE IF NOT USING PROVIDED SAMPLE DATA
par_path=""
ref_fa="dm6.fa" # path to the reference genome .fasta file
declare -a parents=("par1" "par2") # names of parent fly lines

echo $(date +%m-%d-%y)
echo $(date +"%T")
echo

for par in "${parents[@]}"
do
    if [ "${steps:0:1}" = "y" ]
        then
        # ALIGN FILES
        echo -----------STARTING $par ALIGNMENT------------
        bwa mem -R "@RG\tID:$par\tSM:$par" "$ref_fa" "$par_path""$par"_R1.fastq "$par_path""$par"_R2.fastq > "$par"_"$ref".sam
        echo -----------CONVERTING FROM SAM TO BAM------------
        samtools view -S -h -b "$par"_"$ref".sam > "$par"_"$ref".bam
        echo
    fi

    if [ "${steps:1:1}" = "y" ]
        then
        # SORT BAM FILE
        echo -----------SORTING BAM FILES: "$par"------------
        samtools sort "$par"_"$ref".bam > "$par"_"$ref"_first_sort.bam
    fi

    if [ "${steps:2:1}" = "y" ]
        then
        # EXTRACT READS LOCATED IN MICROSATELLITES FROM BAM FILE
        echo -----------EXTRACTING BAM FILES $par------------
        python3 extract_microsat_regions_from_bam_file.py -b "$par"_"$ref"_first_sort.bam -l "$par" 
    fi

    if [ "${steps:3:1}" = "y" ]
        then
        # SORT AND INDEX BAM FILE
        echo -----------SORTING EXTRACTED BAM FILES $par------------
        samtools sort "$par"_"$ref"_extracted.bam > "$par"_"$ref"_extracted.sort.bam
        echo -----------INDEXING $par------------
        samtools index "$par"_"$ref"_extracted.sort.bam
    fi

    if [ "${steps:4:1}" = "y" ]
        then
        # CREATE COVERAGE FILES
        echo -----------STARTING COVERAGE $par------------
        bedtools genomecov -bg -ibam "$par"_"$ref"_extracted.sort.bam > "$par"_"$ref"_coverage.txt
    fi
done

echo -----------DONE WITH COVERAGE $par------------