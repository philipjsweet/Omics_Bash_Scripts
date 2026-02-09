## Goal: Run the GATK 4 Microbes mutation calling tool
## Inputs: XX_s_h_rd.bam resulting from GATK_prep.sh that has a proper header and PCR/Seq duplicates removed
## Reference Genome: Species specific .fna

## Steps
## 1. Index Ref. Genome with samtools
## 2. Index Ref. Genome with Picard
## 3. Index input .bam with samtools
## 4. Run GATK mutect2 
## 5. Use microbes specific GATK filtering
## 6. Move files back to AWS

## Ec2 Chip: m5.4xlarge w/16 vCPUs and 200Gb
## Ec2 Chip: m1.xlarge w/4 vCPUs, 34gb RAM and 200GB attached (not needed)
## EC2 CHi[: m5.xlarge 3/4 vCPU and 16gb of RAM and 200GB attached
## Modules needed
# conda create --name mutect2_env bioconda::picard bioconda::gatk4 bioconda::samtools


## Variables
REF_PATH='s3://bn-seq-data/Ref_Genomes/ScY55_1034/GCA_903819135.2_Y55_Assembly_v2.1_genomic_1034_7.fna'

RUNNER_FILE='runner_GATK_PP_scy55.sh'

SOURCE='s3://bn-seq-data/bn1_pp/wgs/scy55/GATK_Inputs/'

OUTPUT_FOLDER='s3://bn-seq-data/bn1_pp/wgs/scy55/GATK_Outputs/'

REF_NAME='GCA_903819135.2_Y55_Assembly_v2.1_genomic_1034_7.fna'

FOLDER="gatk_scY55"

## Code
printf "" > $RUNNER_FILE

aws s3 ls $SOURCE > test.txt # get files in folder

awk '$4 ~ /\_s_h_rd.bam$/ {print $4}' test.txt > files.txt # UPDATE FOR FILE NAMES

echo "mkdir reads" >> $RUNNER_FILE
echo "mkdir ref_genome" >> $RUNNER_FILE
echo "mkdir $FOLDER" >> $RUNNER_FILE

echo "aws s3 cp $REF_PATH /home/ec2-user/ref_genome" >> $RUNNER_FILE
echo "samtools faidx ref_genome/$REF_NAME" >> $RUNNER_FILE
echo "picard CreateSequenceDictionary -R ref_genome/$REF_NAME"  >> $RUNNER_FILE

while read file; do

RAW_VCF=${file/.bam/.vcf}
FILTERED_VCF=${RAW_VCF/.vcf/_filt.vcf}
BAMOUT=${file/.bam/_gatk.bam}

## Header
echo "## $file " >> $RUNNER_FILE

## Get Files
echo "aws s3 cp $SOURCE$file /home/ec2-user/reads" >> $RUNNER_FILE
echo "samtools index reads/$file" >> $RUNNER_FILE

#### Run Mutect2 ######
echo "gatk Mutect2 -R ref_genome/$REF_NAME -I reads/$file -O $FOLDER/$RAW_VCF -bamout $FOLDER/$BAMOUT" >> $RUNNER_FILE

## Default Filters
echo "gatk FilterMutectCalls --microbial-mode -R ref_genome/$REF_NAME -V $FOLDER/$RAW_VCF -O $FOLDER/$FILTERED_VCF" >> $RUNNER_FILE

## Results back to S3

echo "aws s3 mv $FOLDER/ $OUTPUT_FOLDER --recursive" >> $RUNNER_FILE

## Clear files
echo "rm reads/$file" >> $RUNNER_FILE

done < files.txt


