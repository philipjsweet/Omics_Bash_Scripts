## GOAL: Run BWA 
## 1. Get Reference Genome 
## 2. Index Genome 
## 3. Get Reads
## 4. Align Reads with BWA
## 5. Sort the Reads (reduce size)
## 6. Get FlagStats
## 7. Transfer .bam to S3
## 8. At end, zip all Flagstats and transfer to S3

## Moduels Needed
## 1. BWA
## 2. SAMTOOLS
## 3. multiQC

## Variables 
READ_SOURCE='s3://bn-seq-data/bn2_pb/ScB_2022/trimmed/'
REF_GENOME='s3://bn-seq-data/Ref_Genomes/ScB_2022/'
RUNNER_FILE='bwa_runner_pb_ScY55_1034.sh'
INPUT_FOLDER='/home/ec2-user/reads/'
OUTPUT_FOLDER='s3://bn-seq-data/bn2_pb/ScB_2022/aligned/'
THREADS='16'

printf '' > $RUNNER_FILE

## Set Up Ref Genome Folder
mkdir Ref_Genome

aws s3 cp $REF_GENOME  /home/ec2-user/Ref_Genome/ --recursive

bwa index  /home/ec2-user/Ref_Genome/*.fna

## Get Trimmed Reads 
aws s3 ls $READ_SOURCE > test.txt # get files in folder

awk '$4 ~ /\R1_raw_val_1.fq.gz$/ {print $4}' test.txt > files.txt # UPDATE FOR FILE NAMES

echo "mkdir reads" >> $RUNNER_FILE
echo "mkdir flag_logs" >> $RUNNER_FILE
## Running BWA

while read R1; do

echo "# $R1"  >> $RUNNER_FILE

temp1=$(echo $R1 | awk -F "R1" '{print $1 "R2" $2}')
R2=$(echo $temp1 | awk -F "val_1" '{print $1 "val_2" $2}')

temp2=${R1/R1_raw_val_1/_bwaMem} # replace prefix

OUTPUT=${temp2/.fq.gz/.bam}  # replace file type

SORTED=${OUTPUT/.bam/_sorted.bam}  # replace file type

LOG=${OUTPUT/.bam/_log.tsv} ## Flagstat Log


# Move file from S3 to EC2
echo "aws s3 cp $READ_SOURCE$R1 /home/ec2-user/reads" >> $RUNNER_FILE
echo "aws s3 cp $READ_SOURCE$R2 /home/ec2-user/reads" >> $RUNNER_FILE

# Align RF
echo "bwa mem -t $THREADS /home/ec2-user/Ref_Genome/*.fna /home/ec2-user/reads/$R1 /home/ec2-user/reads/$R2 > reads/$OUTPUT" >> $RUNNER_FILE

echo "samtools sort --threads $THREADS reads/$OUTPUT -o reads/$SORTED" >> $RUNNER_FILE
echo "samtools flagstat --threads $THREADS reads/$SORTED -O tsv > flag_logs/$LOG" >> $RUNNER_FILE
echo "samtools index --threads $THREADS reads/$SORTED" >> $RUNNER_FILE
# Move back to S3

echo "rm /home/ec2-user/reads/*fq.gz" >> $RUNNER_FILE
echo "rm /home/ec2-user/reads/$OUTPUT" >> $RUNNER_FILE
echo "aws s3 mv /home/ec2-user/reads/ $OUTPUT_FOLDER --recursive" >> $RUNNER_FILE

done < files.txt

echo "zip -r flag_logs.zip flag_logs" >> $RUNNER_FILE
echo "aws s3 mv /home/ec2-user/flag_logs.zip $OUTPUT_FOLDER" >> $RUNNER_FILE

#echo "rm test.txt" >> $RUNNER_FILE

## Only got up to ~10% of memory. Started at 10:06, finsihed at 15! so ~10min per file. So fast. 98.94% mapped








