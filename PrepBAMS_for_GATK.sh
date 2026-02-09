## Goal: Edit the bam header and remove PCR duplicates and record changes
## https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
## ID = Read Group ID: Used when dealing with multiplexing. Single library preparation derived from a single biological sample
##       was run on a single lane of a flow cell, all the reads from that lane run belong to the same read group
## PU = Platform Unit: The PU holds three types of information, the {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}. 
##       The {FLOWCELL_BARCODE} refers to the unique identifier for a  particular flow cell. 
##       The {LANE} indicates the lane of the flow cell 
##       The {SAMPLE_BARCODE} is a sample/library-specific identifier. 
##       Although the PU is not required by GATK but takes precedence over ID for base recalibration if it is present.
##       In the example shown earlier, two read group fields, ID and PU, appropriately differentiate flow cell lane, marked by .2, a factor that contributes to batch effects.

## SM = Sample: The name of the sample sequenced in this read group. 
##      GATK tools treat all read groups with the same SM value as containing sequencing data for the same sample, 
##      and this is also the name that will be used for the sample column in the VCF file.

## PL = Platform/technology used to produce the read: Valid values: ILLUMINA, SOLID, LS454, HELICOS and PACBIO.

## LB = DNA preparation library identifier: MarkDuplicates uses the LB field to determine which read groups might contain molecular duplicates

## Modules: conda create -n gatk_prep_env bioconda::samtools bioconda::picard
# Path
SOURCE='s3://bn-seq-data/bn2_pb/KL_04F/aligned/'
OUTPUT='s3://bn-seq-data/bn2_pb/KL_04F/GATK_Inputs/'
## Runner Name
RUNNER_FILE='GATK_prep_runner_PB_KL.sh'
## Library Prep
LIB='PB'

printf "" > $RUNNER_FILE

## Get Trimmed Reads 
aws s3 ls $SOURCE > test.txt # get files in folder

awk '$4 ~ /\__bwaMem_sorted.bam$/ {print $4}' test.txt > files.txt # UPDATE FOR FILE NAMES

echo "mkdir reads" >> $RUNNER_FILE
echo "mkdir picard_logs" >> $RUNNER_FILE

while read BAMS; do

SM=$(echo $BAMS | awk -F "_bwaMem_sorted.bam" '{print $1 "BWA" $2}') 
OUTPUT_1=${BAMS/sorted.bam/s_h.bam} # replace prefix
OUTPUT_2=${BAMS/sorted.bam/s_h_rd.bam}

## Header
echo "## $SM " >> $RUNNER_FILE

## Get Files
echo "aws s3 cp $SOURCE$BAMS /home/ec2-user/reads" >> $RUNNER_FILE
## Set Headers
echo "picard AddOrReplaceReadGroups -I reads/$BAMS -O \
 reads/$OUTPUT_1 \
 -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM $SM" >> $RUNNER_FILE

## Remove Duplicates
echo "picard MarkDuplicates -I reads/$OUTPUT_1 \
    -O reads/$OUTPUT_2 \
    -M picard_logs/marked_dup_metrics.$SM.txt --REMOVE_DUPLICATES true" >> $RUNNER_FILE

## Delete Input
echo "rm reads/$BAMS" >> $RUNNER_FILE
## Move back to S3 
echo "aws s3 mv reads $OUTPUT --recursive" >> $RUNNER_FILE

done < files.txt

echo "zip -r picard_logs.zip picard_logs" >> $RUNNER_FILE
echo "aws s3 mv /home/ec2-user/picard_logs.zip $OUTPUT" >> $RUNNER_FILE

## Takes ~10min per file on a c6g.large. Can't thread the picard steps

