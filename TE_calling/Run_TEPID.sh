#Run Trimmomnatic
# for example on the Ko3_1 files:
trimmomatic PE Ko3_1.fq.gz Ko3_2.fq.gz Ko3_R1_paired.fastq.gz Ko3_R1_unpaired.fastq.gz Ko3_R2_paired.fastq.gz Ko3_R2_unpaired.fastq.gz ILLUMINACLIP:/cluster/work/gdc/people/rhorvath/Scripts/Adapter_Gordon_2017.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Run FastQC
# for example on the Ko3_R1_paired.fastq.gz file:
fastqc Ko3_R1_paired.fastq.gz

# Run TEPID
# for example on the Ko3 file:
tepid-map -x /cluster/scratch/rhorvath/reference_genome/bowtie2_ref_index/Bdistachyon_556_v3.0 -y /cluster/scratch/rhorvath/reference_genome/yaha_ref_index/Bdistachyon_556_v3.0.X15_01_65525S -p 8 -s 300 -n Ko3 -1 /cluster/scratch/rhorvath/FastQ_files/Ko3_R1_paired.fastq.gz -2 /cluster/scratch/rhorvath/FastQ_files/Ko3_R2_paired.fastq.gz -z

