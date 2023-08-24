#Run Trimmomnatic
# for example on the Ko3_1 files:
trimmomatic PE Ko3_1.fq.gz Ko3_2.fq.gz Ko3_R1_paired.fastq.gz Ko3_R1_unpaired.fastq.gz Ko3_R2_paired.fastq.gz Ko3_R2_unpaired.fastq.gz ILLUMINACLIP:/cluster/work/gdc/people/rhorvath/Scripts/Adapter_Gordon_2017.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# Merge paired and unpaired files:
zcat Ko3_R1_paired.fastq.gz Ko3_R1_unpaired.fastq.gz > Ko3_R1.fastq.gz
zcat Ko3_R2_paired.fastq.gz Ko3_R2_unpaired.fastq.gz > Ko3_R2.fastq.gz

# Run FastQC
# for example on the Ko3_R1.fastq.gz file:
fastqc Ko3_R1.fastq.gz

# Run TEPID
# for example on the Ko3 file:
tepid-map -x /cluster/scratch/rhorvath/reference_genome/bowtie2_ref_index/Bdistachyon_556_v3.0 -y /cluster/scratch/rhorvath/reference_genome/yaha_ref_index/Bdistachyon_556_v3.0.X15_01_65525S -p 8 -s 300 -n Ko3 -1 /cluster/scratch/rhorvath/FastQ_files/Ko3_R1.fastq.gz -2 /cluster/scratch/rhorvath/FastQ_files/Ko3_R2.fastq.gz -z
tepid-discover -p 1 -n v -c Ko3.bam -s Ko3.split.bam -t /cluster/scratch/rhorvath/TE_annotation/Bdistachyon_314_v3.0.fa.mod.EDTA.TEanno.bed

# refine step after all samples were run
merge_insertions.py -f insertions 
merge_deletions.py -f deletions

tepid-refine -i /cluster/work/gdc/people/rhorvath/TEPID_output/insertions.bed -d /cluster/work/gdc/people/rhorvath/TEPID_output/deletions.bed -p 4 -t /cluster/work/gdc/people/rhorvath/Reference/TE_annotation/Bdistachyon_314_v3.0.fa.mod.EDTA.TEanno.bed -n Ko3 -c Ko3.bam -s Ko3.split.bam -a /cluster/work/gdc/people/rhorvath/TEPID_output/ALL_Sample_names.txt




