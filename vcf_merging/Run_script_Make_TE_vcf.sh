## run the Make_TE_vcf.R script:

bsub -n 1 -W 24:00 "Rscript --vanilla Make_clean_TE_vcf_file.R 1 10000"
bsub -n 1 -W 24:00 "Rscript --vanilla Make_clean_TE_vcf_file.R 10001 20000"
bsub -n 1 -W 24:00 "Rscript --vanilla Make_clean_TE_vcf_file.R 20001 30000"
bsub -n 1 -W 24:00 "Rscript --vanilla Make_clean_TE_vcf_file.R 30001 40000"
bsub -n 1 -W 24:00 "Rscript --vanilla Make_clean_TE_vcf_file.R 40001 50000"
bsub -n 1 -W 24:00 "Rscript --vanilla Make_clean_TE_vcf_file.R 50001 60000"
bsub -n 1 -W 24:00 "Rscript --vanilla Make_clean_TE_vcf_file.R 60001 70000"
bsub -n 1 -W 24:00 "Rscript --vanilla Make_clean_TE_vcf_file.R 70001 80000"
bsub -n 1 -W 24:00 "Rscript --vanilla Make_clean_TE_vcf_file.R 80001 90000"
bsub -n 1 -W 24:00 "Rscript --vanilla Make_clean_TE_vcf_file.R 90001 100000"
bsub -n 1 -W 24:00 "Rscript --vanilla Make_clean_TE_vcf_file.R 100001 110000"
bsub -n 1 -W 24:00 "Rscript --vanilla Make_clean_TE_vcf_file.R 110001 120000"
bsub -n 1 -W 24:00 "Rscript --vanilla Make_clean_TE_vcf_file.R 120001 128972"

gunzip Full_clean_TE_call_TEPID_1_10000.vcf.gz
gunzip Full_clean_TE_call_TEPID_10001_20000.vcf.gz
gunzip Full_clean_TE_call_TEPID_20001_30000.vcf.gz
gunzip Full_clean_TE_call_TEPID_30001_40000.vcf.gz
gunzip Full_clean_TE_call_TEPID_40001_50000.vcf.gz
gunzip Full_clean_TE_call_TEPID_50001_60000.vcf.gz
gunzip Full_clean_TE_call_TEPID_60001_70000.vcf.gz
gunzip Full_clean_TE_call_TEPID_70001_80000.vcf.gz
gunzip Full_clean_TE_call_TEPID_80001_90000.vcf.gz
gunzip Full_clean_TE_call_TEPID_90001_100000.vcf.gz
gunzip Full_clean_TE_call_TEPID_100001_110000.vcf.gz
gunzip Full_clean_TE_call_TEPID_110001_120000.vcf.gz
gunzip Full_clean_TE_call_TEPID_120001_128972.vcf.gz

head -14 Full_clean_TE_call_TEPID_1_10000.vcf > header_TE_vcf
tail -10000 Full_clean_TE_call_TEPID_1_10000.vcf > only_call_Full_clean_TE_call_TEPID_1_10000.vcf
tail -10000 Full_clean_TE_call_TEPID_10001_20000.vcf > only_call_Full_clean_TE_call_TEPID_10001_20000.vcf
tail -10000 Full_clean_TE_call_TEPID_20001_30000.vcf > only_call_Full_clean_TE_call_TEPID_20001_30000.vcf
tail -10000 Full_clean_TE_call_TEPID_30001_40000.vcf > only_call_Full_clean_TE_call_TEPID_30001_40000.vcf
tail -10000 Full_clean_TE_call_TEPID_40001_50000.vcf > only_call_Full_clean_TE_call_TEPID_40001_50000.vcf
tail -10000 Full_clean_TE_call_TEPID_50001_60000.vcf > only_call_Full_clean_TE_call_TEPID_50001_60000.vcf
tail -10000 Full_clean_TE_call_TEPID_60001_70000.vcf > only_call_Full_clean_TE_call_TEPID_60001_70000.vcf
tail -10000 Full_clean_TE_call_TEPID_70001_80000.vcf > only_call_Full_clean_TE_call_TEPID_70001_80000.vcf
tail -10000 Full_clean_TE_call_TEPID_80001_90000.vcf > only_call_Full_clean_TE_call_TEPID_80001_90000.vcf
tail -10000 Full_clean_TE_call_TEPID_90001_100000.vcf > only_call_Full_clean_TE_call_TEPID_90001_100000.vcf
tail -10000 Full_clean_TE_call_TEPID_100001_110000.vcf > only_call_Full_clean_TE_call_TEPID_100001_110000.vcf
tail -10000 Full_clean_TE_call_TEPID_110001_120000.vcf > only_call_Full_clean_TE_call_TEPID_110001_120000.vcf
tail -10000 Full_clean_TE_call_TEPID_120001_128972.vcf | head -8972 > only_call_Full_clean_TE_call_TEPID_120001_128972.vcf

## merge the vcf files:
cat header_TE_vcf only_call_Full_clean_TE_call_TEPID_1_10000.vcf only_call_Full_clean_TE_call_TEPID_100001_110000.vcf only_call_Full_clean_TE_call_TEPID_40001_50000.vcf only_call_Full_clean_TE_call_TEPID_10001_20000.vcf only_call_Full_clean_TE_call_TEPID_50001_60000.vcf only_call_Full_clean_TE_call_TEPID_120001_128972.vcf only_call_Full_clean_TE_call_TEPID_60001_70000.vcf only_call_Full_clean_TE_call_TEPID_120001_130000.vcf only_call_Full_clean_TE_call_TEPID_70001_80000.vcf only_call_Full_clean_TE_call_TEPID_80001_90000.vcf only_call_Full_clean_TE_call_TEPID_20001_30000.vcf only_call_Full_clean_TE_call_TEPID_90001_100000.vcf only_call_Full_clean_TE_call_TEPID_30001_40000.vcf > Full_clean_TE_call_TEPID_unsorted.vcf

## sort vcf:
vcf-sort Full_clean_TE_call_TEPID_unsorted.vcf > Full_clean_TE_call_TEPID.vcf

## correct some formatin issues:
sed -i -e 's/,T\t/\t/g' Full_clean_TE_call_TEPID.vcf

## Select the same samples as in the SNP vcf:
gatk IndexFeatureFile -F ./Full_clean_TE_call_TEPID.vcf
gatk SelectVariants -R ./Bdistachyon_556_v3.0.fa -V ./Full_clean_TE_call_TEPID.vcf --exclude-sample-expressions SRR1800522 --exclude-sample-expressions Bd1-1 --exclude-sample-expressions Bd21-3 --exclude-sample-expressions BdTR8i --exclude-sample-expressions BdTR3C --exclude-sample-expressions Bd18-1 --exclude-sample-expressions Cef2 --exclude-sample-expressions Pob1 --exclude-sample-expressions D26 --exclude-sample-expressions ABR8 --exclude-sample-expressions Bd21Control_1 --exclude-sample-expressions Bd21Control_2 --exclude-sample-expressions Bd21Control_3 --exclude-sample-expressions Bd21_Skalska --exclude-sample-expressions Bd21_Stritt --exclude-sample-expressions D73 --exclude-non-variants true -O Bdis326_clean_TE_call_TEPID.vcf
