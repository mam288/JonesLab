#!/bin/bash

# Input Parameters
steps=$1 # which processing steps to perform Align: Sort: Extract: Create Coverage Files: Call Variants: Split Variant Files: Process Microsatellite Files
trunc=$2 #how many lines to truncate the tables to
save=$3 

# Set variables
# CHANGE THESE IF DESIRED
ct=5 #coverage threshold
ad=2 #allele depth threshold
mq=0 #mapping quality threshold
b=5 # length of prefix and suffix sequences to pull
m=3 # factor to multiply the buffer by when getting the prefix and suffix
sd=1500 # how many bp apart the start distances can be to match the rows
covtype=bg # type of coverage file - change to d if calling bedtools genome_cov with -d options instead of -bg
par1_name="par1" # name of the first parent fly line
par2_name="par2" # name of the second parent line line
offspring_name="offspring" # name of the offspring fly line
ref="dm6" # ref genome name

# Set file paths 
# CHANGE THESE IF NOT USING PROVIDED SAMPLE DATA
fastq_1="offspring_R1.fastq" # path to first fastq file
fastq_2="offspring_R2.fastq" # path to second fastq file
p1_microsats="par1_microsats.td" # microsatellite file for parent 1
p2_microsats="par2_microsats.td" # microsatellite file for parent 2
p1_consensus="par1_consensus.fa" # consensus file for parent 1
p2_consensus="par2_consensus.fa" # consensus file for parent 2
ref_fa="dm6.fa" # path to the reference genome .fasta file
py_script="../scripts/microsatellites_compare_parent_offspring.py" # python script to analyze microsatellites

echo $(date +%m-%d-%y)
echo $(date +"%T")
echo
echo ------Input Parameters: "$steps" "$trunc" "$save"------
echo 

# ALIGN FILES
if [ "${steps:0:1}" = "y" ]
then
    echo -----------STARTING ALIGNMENT: "$offspring_name"------------

    bwa mem -R "@RG\tID:$offspring_name\tSM:$offspring_name" "$ref_fa" "$fastq_1" "$fastq_2" \
        > "$offspring_name"_"$ref".sam # align files and add the name of the offspring line as the read group
    echo -----------FINISHED ALIGNMENT: "$fly_line"------------
    samtools view -S -h -b "$offspring_name"_"$ref".sam > "$offspring_name"_"$ref".bam
    echo -----------FINISHED SAM TO BAM CONVERSION: "$offspring_name"------------
fi
echo

# SORT BAM FILE
if [ "${steps:1:1}" = "y" ]
then
    #SORT BAM FILE
    echo -----------SORTING BAM FILES: "$offspring_name"------------
    echo samtools sort "$offspring_name"_"$ref".bam > "$offspring_name"_"$ref"_first_sort.bam 
    samtools sort "$offspring_name"_"$ref".bam > "$offspring_name"_"$ref"_first_sort.bam 

    echo -----------FINISHED FIRST SORTING: "$offspring_name"------------
fi

# EXTRACT READS LOCATED IN MICROSATELLITES FROM BAM FILE, SORT AND INDEX BAM FILES
if [ "${steps:2:1}" = "y" ]
then
    # EXTRACT READS LOCATED IN MICROSATELLITES FROM BAM FILE
    echo -----------EXTRACTING BAM FILES: "$offspring_name"------------

    python3 extract_microsat_regions_from_bam_file.py \
        -b "$offspring_name"_"$ref"_first_sort.bam -l "$offspring_name" 

    echo -----------FINISHED EXTRACTION: "$offspring_name"------------

    # SORT AND INDEX BAM FILE
    echo -----------SORTING EXTRACTED BAM FILES: "$offspring_name"------------
    samtools sort "$offspring_name"_"$ref"_extracted.bam > "$offspring_name"_"$ref"_extracted.sort.bam
    echo -----------FINISHED SORTING: "$offspring_name"------------

    echo -----------INDEXING: "$offspring_name"------------
    samtools index "$offspring_name"_"$ref"_extracted.sort.bam
    echo -----------FINISHED INDEXING: "$offspring_name"------------
fi
echo  

# CREATE COVERAGE FILES
if [ "${steps:3:1}" = "y" ]
then
    echo -----------STARTING COVERAGE: "$fly_line_name"------------
    bedtools genomecov -bg -ibam "$offspring_name"_"$ref"_extracted.sort.bam > "$offspring_name"_"$ref"_coverage.txt
    echo -----------DONE WITH COVERAGE "$fly_line_name"------------
fi

# CALL VARIANTS WITH SELECTED VARIANT CALLER
if [ "${steps:4:1}" = "y" ]
then
    echo -----------STARTING VARIANT CALLING: "$fly_line_name" "$vc"------------

    freebayes -f "$ref_fa" --min-alternate-count 1 --min-alternate-fraction 0 -b "$offspring_name"_"$ref"_extracted.sort.bam \
        "$par1_name"_"$ref"_extracted.sort.bam "$par2_name"_"$ref"_extracted.sort.bam  > "$offspring_name"_with_parents.vcf
	echo ------------DONE WITH VARIANT CALLING------------
fi

# SPLIT VARIANT FILE
if [ "${steps:5:1}" = "y" ]
then
    file_name="$offspring_name"_with_parents.vcf
    rm $file_name.gz # remove compressed file if it is there, so a new one can be made
    bgzip $file_name # compress file
    echo -----------INDEXING MASTER VCF FILE: "$fly_line_name"------------
    bcftools tabix -p vcf -f $file_name.gz
    echo -----------DONE INDEXING: "$fly_line_name"------------
    
    echo -----------SEPARATING VCF FILE: "$fly_line_name"------------
    echo 
    samples=$(bcftools query -l $file_name.gz)
    echo samples $samples
    all_files=0
    for sample in `bcftools query -l $file_name.gz`
    do
        echo
        echo -----------SEPARATING SAMPLE: $sample------------
        new_file_name=${file_name/.vcf/.$sample.vcf.gz}
        echo new file $new_file_name
        bcftools view -c1 -Oz -s $sample -o $new_file_name $file_name.gz
        echo -----------STARTING TABIX: "$fly_line_name" $sample------------
        bcftools tabix -p vcf -f $new_file_name
        echo -----------FINISHED TABIX: "$fly_line_name" $sample------------

    done
fi


echo  

# PROCESS MICROSATELLITE FILES
if [ "${steps:6:1}" = "y" ]
then
	echo ------------PROCESSING "$fly_line_name"-----------------

    echo Input Parameters: -s "$trunc" -b "$b" -m "$m" -sd "$sd" -ct "$ct" -ad "$ad" -mq "$mq" -covtype "$covtype" -gf "$ref_fa"

    python3 "$py_script" \
    -s "$trunc" -b "$b" -m "$m" -sd "$sd" -ct "$ct" -ad "$ad" -mq "$mq" -covtype "$covtype" -gf "$ref_fa" \
    -gmf "$ref"_microsats.td \
    -vcf "$offspring_name"_with_parents."$par1_name".vcf.gz "$offspring_name"_with_parents."$par2_name".vcf.gz "$offspring_name"_with_parents."$offspring_name".vcf.gz \
    -cf "$par1_name"_"$ref"_coverage.txt  "$par2_name"_"$ref"_coverage.txt "$offspring_name"_"$ref"_coverage.txt \
    -mf "$p1_microsats" "$p2_microsats" \
    -ff "$p1_consensus" "$p2_consensus"  \
    -save "$save"

echo ------------DONE PROCESSING-----------------
fi
echo
echo ------------DONE-----------------

