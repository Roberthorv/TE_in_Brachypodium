## Prepare the genetic map called Bdistachyon_genetic_map.txt
sed -i -e 's/_//g' Bdistachyon_genetic_map.txt

## split genetic map by scaffolds:
grep "chr1\|Chromosome" Bdistachyon_genetic_map.txt > Bdistachyon_genetic_map_chr1.txt
grep "chr2\|Chromosome" Bdistachyon_genetic_map.txt > Bdistachyon_genetic_map_chr2.txt
grep "chr3\|Chromosome" Bdistachyon_genetic_map.txt > Bdistachyon_genetic_map_chr3.txt
grep "chr4\|Chromosome" Bdistachyon_genetic_map.txt > Bdistachyon_genetic_map_chr4.txt
grep "chr5\|Chromosome" Bdistachyon_genetic_map.txt > Bdistachyon_genetic_map_chr5.txt


## remove multiallelic sites from vcf file
gatk SelectVariants -R ./Bdistachyon_556_v3.0.fa -V Bdis332_SNPs_Yann_clean_TE_TEPID.vcf.gz --restrict-alleles-to BIALLELIC -O Bdis326_Biallelic_SNPs_clean_TE_TEPID_GEVA.vcf.gz

## Polerizing the vcf using the C clade:
## select C clade:
gatk SelectVariants -R ./Bdistachyon_556_v3.0.fa -V Bdis326_Biallelic_SNPs_clean_TE_TEPID_GEVA.vcf.gz --sample-expressions D71 --sample-expressions D72 --sample-expressions Cb23 --sample-expressions Cro24 --sample-expressions D65 --sample-expressions D66 --sample-expressions D67 --sample-expressions D68 --sample-expressions D69 --sample-expressions D70 --sample-expressions Mca12 --sample-expressions Mca16 --sample-expressions Mca18 --sample-expressions Mca19 --sample-expressions Mca23 --sample-expressions D25 --sample-expressions Msa11 --sample-expressions Msa16 --sample-expressions Msa6 --sample-expressions Msa9 --sample-expressions San11 --sample-expressions San12 --sample-expressions San25 --sample-expressions Spm16 --sample-expressions SRR4162922 --sample-expressions SRR4162932 --sample-expressions SRR4162933 -O Bdis326_Biallelic_SNPs_clean_TE_TEPID_GEVA_C.vcf.gz

## select only fixed sites in C clade:
zcat Bdis326_Biallelic_SNPs_clean_TE_TEPID_GEVA_C.vcf.gz | grep "AF=0.00\|AF=1.00" | cut -f 1,2 > Bdis326_fixed_sites_GEVA_C.txt
awk '{print $1 "\t" ($2 - 1) "\t" $2}' Bdis326_fixed_sites_GEVA_C.txt > Bdis326_fixed_sites_GEVA_C.bed
zcat Bdis326_Biallelic_SNPs_clean_TE_TEPID_GEVA_C.vcf.gz | grep "AF=1.00" | cut -f 1,2,4,5 > Fixed_alternate_sites_C_for_polarization.txt
awk '{print $1 "\t" ($2 - 1) "\t" $2}' Fixed_alternate_sites_C_for_polarization.txt > Fixed_alternate_sites_C_for_polarization.bed

## select only sites that can be polarized:
zcat Bdis326_Biallelic_SNPs_clean_TE_TEPID_GEVA.vcf.gz | head -61 > head_vcf
bedtools intersect -a Bdis326_Biallelic_SNPs_clean_TE_TEPID_GEVA.vcf.gz -b Bdis326_fixed_sites_GEVA_C.bed -wa > Bdis326_SNPs_TE_for_polarization_body.vcf
cat head_vcf Bdis326_SNPs_TE_for_polarization_body.vcf | bedtools intersect -a - -b Fixed_alternate_sites_C_for_polarization.bed -wa > Fixed_alternate_sites_C_SNPs_TE_body.vcf
cat head_vcf Bdis326_SNPs_TE_for_polarization_body.vcf | bedtools subtract -a - -b Fixed_alternate_sites_C_for_polarization.bed -wa > Fixed_ancestral_sites_C_SNPs_TE_body.vcf

## split vcf for running in paralle
head -50000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | cat head_vcf - > To_fixed_SNPs_TE_part_1.vcf
head -100000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_2.vcf
head -150000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_3.vcf
head -200000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_4.vcf
head -250000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_5.vcf
head -300000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_6.vcf
head -350000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_7.vcf
head -400000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_8.vcf
head -450000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_9.vcf
head -500000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_10.vcf
head -550000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_11.vcf
head -600000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_12.vcf
head -650000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_13.vcf
head -700000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_14.vcf
head -750000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_15.vcf
head -800000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_16.vcf
head -850000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_17.vcf
head -900000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_18.vcf
head -950000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_19.vcf
head -1000000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_20.vcf
head -1050000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_21.vcf
head -1100000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_22.vcf
head -1150000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_23.vcf
head -1200000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_24.vcf
head -1250000 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -50000 | cat head_vcf - > To_fixed_SNPs_TE_part_25.vcf
head -1279695 Fixed_alternate_sites_C_SNPs_TE_body.vcf | tail -29695 | cat head_vcf - > To_fixed_SNPs_TE_part_26.vcf


## polariz the different parts using the Polarize_vcf_file.R script:
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 1"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 2"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 3"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 4"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 5"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 6"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 7"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 8"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 9"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 10"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 11"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 12"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 13"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 14"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 15"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 16"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 17"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 18"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 19"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 20"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 21"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 22"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 23"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 24"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 25"
bsub -n 1 -W 24:00 "Rscript --vanilla ./Polarize_vcf_file.R 26"

zcat Pol_SNPs_TE_part_1.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_1_body.vcf
zcat Pol_SNPs_TE_part_2.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_2_body.vcf
zcat Pol_SNPs_TE_part_3.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_3_body.vcf
zcat Pol_SNPs_TE_part_4.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_4_body.vcf
zcat Pol_SNPs_TE_part_5.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_5_body.vcf
zcat Pol_SNPs_TE_part_6.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_6_body.vcf
zcat Pol_SNPs_TE_part_7.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_7_body.vcf
zcat Pol_SNPs_TE_part_8.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_8_body.vcf
zcat Pol_SNPs_TE_part_9.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_9_body.vcf
zcat Pol_SNPs_TE_part_10.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_10_body.vcf
zcat Pol_SNPs_TE_part_11.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_11_body.vcf
zcat Pol_SNPs_TE_part_12.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_12_body.vcf
zcat Pol_SNPs_TE_part_13.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_13_body.vcf
zcat Pol_SNPs_TE_part_14.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_14_body.vcf
zcat Pol_SNPs_TE_part_15.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_15_body.vcf
zcat Pol_SNPs_TE_part_16.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_16_body.vcf
zcat Pol_SNPs_TE_part_17.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_17_body.vcf
zcat Pol_SNPs_TE_part_18.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_18_body.vcf
zcat Pol_SNPs_TE_part_19.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_19_body.vcf
zcat Pol_SNPs_TE_part_20.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_20_body.vcf
zcat Pol_SNPs_TE_part_21.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_21_body.vcf
zcat Pol_SNPs_TE_part_22.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_22_body.vcf
zcat Pol_SNPs_TE_part_23.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_23_body.vcf
zcat Pol_SNPs_TE_part_24.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_24_body.vcf
zcat Pol_SNPs_TE_part_25.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_25_body.vcf
zcat Pol_SNPs_TE_part_26.vcf.gz | sed '1,61d' > Pol_SNPs_TE_part_26_body.vcf


## merge and sort vcf 
cat head_vcf Fixed_ancestral_sites_C_SNPs_TE_body.vcf Pol_SNPs_TE_part_*_body.vcf | vcf-sort > Pol_SNPs_TE.vcf

## phase the vcf
sed -e 's@1/1@1|1@g' Pol_SNPs_TE.vcf > Phased_Pol_SNPs_TE.vcf
sed -i -e 's@0/0@0|0@g' Phased_Pol_SNPs_TE.vcf
sed -i -e 's@\./\.@.|.@g' Phased_Pol_SNPs_TE.vcf

gzip Phased_Pol_SNPs_TE.vcf

