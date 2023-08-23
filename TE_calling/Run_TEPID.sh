#Run Trimmomnatic
# for example on the Ko3_1 files:
trimmomatic PE Ko3_1.fq.gz Ko3_2.fq.gz Ko3_R1_paired.fastq.gz Ko3_R1_unpaired.fastq.gz Ko3_R2_paired.fastq.gz Ko3_R2_unpaired.fastq.gz ILLUMINACLIP:/cluster/work/gdc/people/rhorvath/Scripts/Adapter_Gordon_2017.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Run FastQC
# for example on the Ko3_R1_paired.fastq.gz file:
fastqc Ko3_R1_paired.fastq.gz

