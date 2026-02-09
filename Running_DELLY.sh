## GOAL:
## 1. Get indexed reference genome
## 2. Get production sample bam file (SP_KP_1)
## 3. Get timpoint bam files
## 4. Run Delly
## 5. Delete indexed reference genome from e2
## 6. Delete bam files from e2
## 7. Move SV calls to S3

## Modules Needed
## 1. Delly
## 2. BCFtools, Matplotlib, Numpy, LaTeX
## 3. Samtools

## Variables
READ_SOURCE='s3://bn-seq-data/bn1_pp/wgs/scb/GATK_Inputs/'
REF_GENOME='/home/ec2-user/reference/GCA_001298375.2_ASM129837v2_genomic_2023_4.fna'
PROD_FILE='production_bams/BG_ScY55_1__bwaMem_s_h_rd.bam'
RUNNER_FILE='delly_allTimepoints_healthyRef_PP_scb_runner.sh'
INPUT_FOLDER='/home/ec2-user/alignments/'
LOCAL_OUTPUT='/home/ec2-user/sv_calls/delly/'
OUTPUT_FOLDER='s3://bn-seq-data/bn1_pp/wgs/scb/Delly_Outputs/'

printf '' > $RUNNER_FILE

# Copy reference genome
aws s3 cp s3://bn-seq-data/Ref_Genomes/ScB_2022/GCA_001298375.2_ASM129837v2_genomic_2023_4.fna /home/ec2-user/reference/ 
samtools faidx $REF_GENOME

# Copy production bam file and create index
aws s3 cp s3://bn-seq-data/bn1_pp/wgs/Desiccation_Study/GATK_Inputs/ScB/BG_ScY55_1__bwaMem_s_h_rd.bam production_bams/
samtools index $PROD_FILE

## Get files
printf '' > files.txt
aws s3 ls $READ_SOURCE > test.txt # get files in folder
awk '$4 ~ /\_s_h_rd.bam$/ {print $4}' test.txt > files.txt # UPDATE FOR FILE NAMES

## Running DELLY
while read BAM; do # Iterate over files.txt to run delly

echo "# $BAM"  >> $RUNNER_FILE # copies bam file name as comment in to runner file

BAM_INDEX=${BAM/.bam/.bam.bai} # replaces bam file type with bai file type
OUTPUT_FILE="${BAM/__bwaMem_s_h_rd.bam/_SVcalls_delly_allTimepoints_healthyRef.bcf}" # Replace prefix to fit SV calls
VCF_OUTPUT="${BAM/__bwaMem_s_h_rd.bam/_SVcalls_delly_allTimepoints_healthyRef.vcf}"

# Move file from S3 to EC2
echo "aws s3 cp $READ_SOURCE$BAM $INPUT_FOLDER" >> $RUNNER_FILE # Copies bam file
echo "samtools index $INPUT_FOLDER$BAM $INPUT_FOLDER$BAM_INDEX" >> $RUNNER_FILE # Creates bam.bai file

# SV Calls with Delly
echo "delly call -o $LOCAL_OUTPUT$OUTPUT_FILE -g $REF_GENOME $INPUT_FOLDER$BAM $PROD_FILE" >> $RUNNER_FILE #Run Delly
## CHECK OUTPUT FILE NAME!

#Make vcf file
echo "bcftools view $LOCAL_OUTPUT$OUTPUT_FILE > $LOCAL_OUTPUT$VCF_OUTPUT" >> $RUNNER_FILE

# Move back to S3
echo "rm /home/ec2-user/alignments/*.bam" >> $RUNNER_FILE # Remove bam file
echo "rm /home/ec2-user/alignments/*.bai" >> $RUNNER_FILE # Remove bai file
echo "aws s3 mv /home/ec2-user/sv_calls/delly $OUTPUT_FOLDER --recursive" >> $RUNNER_FILE # Move output vcf to s3

done < files.txt



