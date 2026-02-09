## GOAL: This script will look in a user defined S3 folder...
## 1. for the  fastq files to be trimmed
## 2. get the path to the files
## 3. One by one copy the files to the running drive
## 4. Make a folder called "fastqc" for sumamry files
## 5. Run the trimmer command 
## 6. Save the summary file
## 7. Transfer the trimmed file back to a user defined S3
## 8. Delete the raw file

## Modules Needed 
## 1. cutadapt
## 2. fastQC
## 3. conda install trim-galore
## 4. multiQC
## All found in conda fastqc_env 

## Variables 
RUNNER_FILE='trimmer_runner_KL_04F_pb.sh'


## Pathways 
INPUTS='s3://bn-seq-data/bn2_pb/KL_04F/raw_reads/'
OUTPUT='s3://bn-seq-data/bn2_pb/KL_04F/trimmed'

## Get files
aws s3 ls $INPUTS> test.txt # get files in folder

awk '$4 ~ /\R1_raw.fastq.gz$/ {print $4}' test.txt > files.txt # select for correct file type (R1)

mkdir reads
mkdir fastqc

printf "" > $RUNNER_FILE

while read R1_FILE; do

## Section Header
echo "# $R1_FILE"  >> $RUNNER_FILE 

## Paired File
R2_FILE=$(echo $R1_FILE | awk -F "R1" '{print $1 "R2" $2}') 

# Move file from S3 to EC2
echo "aws s3 cp $INPUTS$R1_FILE /home/ec2-user/reads" >> $RUNNER_FILE
echo "aws s3 cp $INPUTS$R2_FILE /home/ec2-user/reads" >> $RUNNER_FILE

# trim file
echo "trim_galore --paired --illumina --fastqc --output_dir /home/ec2-user/fastqc  /home/ec2-user/reads/$R1_FILE /home/ec2-user/reads/$R2_FILE" >> $RUNNER_FILE

# Move back to S3
echo "aws s3 mv /home/ec2-user/fastqc $OUTPUT --recursive" >> $RUNNER_FILE

# Delete files from EC2
echo "rm -f reads/$R1_FILE" >> $RUNNER_FILE
echo "rm -f reads/$R2_FILE" >> $RUNNER_FILE

done < files.txt

rm files.txt
