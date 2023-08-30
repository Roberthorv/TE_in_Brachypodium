## preparing the data

## split vcf in chromosomes:
gatk IndexFeatureFile -F Phased_Pol_SNPs_TE.vcf

gatk SelectVariants -R ./Bdistachyon_556_v3.0.fa -V Phased_Pol_SNPs_TE.vcf -L Bd1 -O Phased_Pol_SNPs_TE_Bd1.vcf
gatk SelectVariants -R ./Bdistachyon_556_v3.0.fa -V Phased_Pol_SNPs_TE.vcf -L Bd2 -O Phased_Pol_SNPs_TE_Bd2.vcf
gatk SelectVariants -R ./Bdistachyon_556_v3.0.fa -V Phased_Pol_SNPs_TE.vcf -L Bd3 -O Phased_Pol_SNPs_TE_Bd3.vcf
gatk SelectVariants -R ./Bdistachyon_556_v3.0.fa -V Phased_Pol_SNPs_TE.vcf -L Bd4 -O Phased_Pol_SNPs_TE_Bd4.vcf
gatk SelectVariants -R ./Bdistachyon_556_v3.0.fa -V Phased_Pol_SNPs_TE.vcf -L Bd5 -O Phased_Pol_SNPs_TE_Bd5.vcf

## split data in different clusters 

for i in  $(ls Phased_Pol_SNPs_TE_Bd*.vcf)
do 
NAME=$(echo $i | sed -e 's/.vcf//g')
gatk SelectVariants -R ./Bdistachyon_556_v3.0.fa -V $i --sample-expressions 1a_29_12 --sample-expressions 1a_29_14 --sample-expressions 1a_29_27 --sample-expressions 1a_30_24 --sample-expressions 1a_30_26 --sample-expressions 1a_30_2 --sample-expressions 1a_30_4 --sample-expressions 1a_31_2 --sample-expressions 1a_31_6 --sample-expressions 1a_32_10 --sample-expressions 1a_32_12 --sample-expressions 1a_32_5 --sample-expressions 1c_25_7 --sample-expressions 1c_26_3 --sample-expressions 1c_26_4 --sample-expressions 1c_27_15 --sample-expressions 1c_27_8 --sample-expressions 1c_28_1 --sample-expressions 1c_28_2 --sample-expressions 1c_34_1 --sample-expressions 1c_34_25 --sample-expressions Adi-10 --sample-expressions Adi-12 --sample-expressions Adi-2 --sample-expressions 18-1 --sample-expressions 21-3_r --sample-expressions 21Control_1 --sample-expressions 21Control_2 --sample-expressions 21Control_3 --sample-expressions 2-3 --sample-expressions 3-1 --sample-expressions TR10C --sample-expressions TR11A --sample-expressions TR11G --sample-expressions TR11I --sample-expressions TR12c --sample-expressions TR13a --sample-expressions TR13C --sample-expressions TR1i --sample-expressions TR2B --sample-expressions TR2G --sample-expressions TR3C --sample-expressions TR5I --sample-expressions TR9K --sample-expressions Bis-1 --sample-expressions Gaz-8 --sample-expressions Kah-1 --sample-expressions Kah-5 --sample-expressions Koz-1 --sample-expressions Koz-3 --sample-expressions ABR8 --sample-expressions SRR4162886 --sample-expressions SRR4162895 --sample-expressions SRR4162905 --sample-expressions SRR4162910 --sample-expressions SRR4162914 --sample-expressions SRR4162917 --sample-expressions SRR4162963 --sample-expressions SRR4163268 --sample-expressions SRR4163269 --sample-expressions SRR4163270 --sample-expressions SRR4163272 --sample-expressions SRR4163273 --sample-expressions SRR4163274 --sample-expressions SRR4163275 --sample-expressions SRR4163276 --sample-expressions SRR4163277 --sample-expressions SRR4163278 --sample-expressions SRR4164006 --sample-expressions SRR4164007 --sample-expressions SRR4164010 --sample-expressions SRR4164011 --sample-expressions SRR4164022 --sample-expressions SRR4164026 --sample-expressions SRR4164047 --sample-expressions SRR4235973 --exclude-non-variants true -O ${NAME}_B_East.vcf
gatk SelectVariants -R ./Bdistachyon_556_v3.0.fa -V $i --sample-expressions 1c_25_14 --sample-expressions 1c_25_15 --sample-expressions 1c_35_1 --sample-expressions 1c_35_7 --sample-expressions 2_14_13 --sample-expressions 2_14_25 --sample-expressions 2_1_7 --sample-expressions 2_5_10 --sample-expressions 2_8_1 --sample-expressions 2_8_3 --sample-expressions 2_9_1 --sample-expressions 2_9_9 --sample-expressions 3_38_12 --sample-expressions 3_38_19 --sample-expressions 3_40_4 --sample-expressions 3_41_1 --sample-expressions 3_41_3 --sample-expressions 3_42_11 --sample-expressions 3_42_19 --sample-expressions 3_54_19 --sample-expressions 3_54_3 --sample-expressions 3_55_13 --sample-expressions 3_55_27 --sample-expressions 4_36_22 --sample-expressions 4_45_19 --sample-expressions 4_47_16 --sample-expressions 4_47_25 --sample-expressions 4_49_22 --sample-expressions 4_49_26 --sample-expressions 4_51_22 --sample-expressions 4_52_15 --sample-expressions 4_52_6 --sample-expressions 1-1 --sample-expressions 29-1 --sample-expressions TR7a --sample-expressions TR8i --sample-expressions D19 --sample-expressions D20 --sample-expressions D21 --sample-expressions D22 --sample-expressions D23 --sample-expressions D24 --sample-expressions D60 --sample-expressions D61 --sample-expressions D62 --sample-expressions D63 --sample-expressions D64 --sample-expressions D51 --sample-expressions D52 --sample-expressions D53 --sample-expressions D54 --sample-expressions D55 --sample-expressions D56 --sample-expressions D42 --sample-expressions D43 --sample-expressions D44 --sample-expressions D45 --sample-expressions D46 --sample-expressions D47 --sample-expressions D48 --sample-expressions D49 --sample-expressions D50 --sample-expressions SRR4162889 --sample-expressions SRR4162891 --sample-expressions SRR4162898 --sample-expressions SRR4162906 --sample-expressions SRR4162908 --sample-expressions SRR4162911 --sample-expressions SRR4162913 --sample-expressions SRR4163271 --sample-expressions SRR4235993 --sample-expressions D57 --sample-expressions D58 --sample-expressions D59 --sample-expressions Tek-2 --sample-expressions Tek-4 --sample-expressions D36 --sample-expressions D37 --sample-expressions D38 --sample-expressions D39 --sample-expressions D40 --sample-expressions D41 --sample-expressions D27 --sample-expressions D28 --sample-expressions D29 --sample-expressions D30 --sample-expressions D31 --sample-expressions D32 --sample-expressions D33 --sample-expressions D34 --sample-expressions D35 --sample-expressions D16 --sample-expressions D17 --sample-expressions D18 --exclude-non-variants true -O ${NAME}_A_East.vcf
gatk SelectVariants -R ./Bdistachyon_556_v3.0.fa -V $i --sample-expressions 2_14_15 --sample-expressions 2_14_20 --sample-expressions 2_20_16 --sample-expressions ABR2 --sample-expressions ABR3 --sample-expressions ABR4 --sample-expressions ABR5 --sample-expressions ABR6 --sample-expressions ABR7 --sample-expressions Arn1 --sample-expressions 30-1 --sample-expressions Foz1 --sample-expressions Jer1 --sample-expressions Lam13 --sample-expressions Lam6 --sample-expressions Luc1 --sample-expressions Mig3 --sample-expressions Mon3 --sample-expressions D10 --sample-expressions D11 --sample-expressions D12 --sample-expressions D13 --sample-expressions D14 --sample-expressions D15 --sample-expressions D1 --sample-expressions D2 --sample-expressions D3 --sample-expressions D4 --sample-expressions D6 --sample-expressions D7 --sample-expressions D8 --sample-expressions D9 --sample-expressions Mur1 --sample-expressions Per1 --sample-expressions RON2 --sample-expressions S8iiC --sample-expressions Shf16 --sample-expressions Shf23 --sample-expressions Shf28 --sample-expressions Shf8 --sample-expressions Sig2 --sample-expressions SRR4162935 --sample-expressions SRR4162938 --sample-expressions SRR4162939 --sample-expressions SRR4162940 --sample-expressions SRR4162941 --sample-expressions SRR4162942 --sample-expressions SRR4162943 --sample-expressions SRR4162944 --sample-expressions SRR4162945 --sample-expressions SRR4162946 --sample-expressions SRR4162947 --sample-expressions SRR4162948 --sample-expressions SRR4162949 --sample-expressions SRR4162950 --sample-expressions SRR4162951 --sample-expressions SRR4162953 --sample-expressions SRR4162954 --sample-expressions SRR4162955 --sample-expressions SRR4162956 --sample-expressions SRR4162959 --sample-expressions SRR4162960 --sample-expressions SRR4162961 --sample-expressions SRR4162962 --sample-expressions SRR4235987 --sample-expressions Tso17 --sample-expressions Tso18 --sample-expressions Tso22 --sample-expressions Tso3 --sample-expressions Tso8 --sample-expressions Uni2 --sample-expressions Vif1 --exclude-sample-expressions D20 --exclude-sample-expressions D21 --exclude-sample-expressions D22 --exclude-sample-expressions D23 --exclude-sample-expressions D24 --exclude-sample-expressions D25 --exclude-sample-expressions D26 --exclude-sample-expressions D27 --exclude-sample-expressions D28 --exclude-sample-expressions D29 --exclude-sample-expressions D30 --exclude-sample-expressions D31 --exclude-sample-expressions D32 --exclude-sample-expressions D33 --exclude-sample-expressions D34 --exclude-sample-expressions D35 --exclude-sample-expressions D36 --exclude-sample-expressions D37 --exclude-sample-expressions D38 --exclude-sample-expressions D39 --exclude-sample-expressions D40 --exclude-sample-expressions D41 --exclude-sample-expressions D42 --exclude-sample-expressions D43 --exclude-sample-expressions D44 --exclude-sample-expressions D45 --exclude-sample-expressions D46 --exclude-sample-expressions D47 --exclude-sample-expressions D48 --exclude-sample-expressions D49 --exclude-sample-expressions D60 --exclude-sample-expressions D61 --exclude-sample-expressions D62 --exclude-sample-expressions D63 --exclude-sample-expressions D64 --exclude-sample-expressions D65 --exclude-sample-expressions D66 --exclude-sample-expressions D67 --exclude-sample-expressions D68 --exclude-sample-expressions D69 --exclude-sample-expressions D70 --exclude-sample-expressions D71 --exclude-sample-expressions D72 --exclude-sample-expressions D16 --exclude-sample-expressions D17 --exclude-sample-expressions D18 --exclude-sample-expressions D19 --exclude-non-variants true -O ${NAME}_B_West.vcf
gatk SelectVariants -R ./Bdistachyon_556_v3.0.fa -V $i --sample-expressions 2_5_21 --sample-expressions ABR9 --sample-expressions Chi2 --sample-expressions Chi3 --sample-expressions Cm10 --sample-expressions Cm12 --sample-expressions Cm13 --sample-expressions Cm16 --sample-expressions Cm18 --sample-expressions Cm2 --sample-expressions Cm3 --sample-expressions Cm4 --sample-expressions Cm7 --sample-expressions Ges19 --sample-expressions Ges1 --sample-expressions Ges26 --sample-expressions Ges30 --sample-expressions Ges5 --sample-expressions Lb11 --sample-expressions Lb13 --sample-expressions Lb16 --sample-expressions Lb18 --sample-expressions Lb1 --sample-expressions Lb21 --sample-expressions Lb23 --sample-expressions Lb3 --sample-expressions Lb7 --sample-expressions Lb8 --sample-expressions Leo2 --sample-expressions Leo7 --sample-expressions Leo8 --sample-expressions Mca10 --sample-expressions Mca26 --sample-expressions Mcp2 --sample-expressions Mcp7 --sample-expressions D5 --sample-expressions Msa27 --sample-expressions Ram6 --sample-expressions Ren14 --sample-expressions Ren19 --sample-expressions Ren22 --sample-expressions Ren4 --sample-expressions Ren5 --sample-expressions San24 --sample-expressions San26 --sample-expressions Sap31 --sample-expressions Sap40 --sample-expressions Sap44 --sample-expressions Sap47 --sample-expressions Sap49 --sample-expressions Sap57 --sample-expressions Sap58 --sample-expressions Sap59 --sample-expressions Sap61 --sample-expressions Sap64 --sample-expressions Sap66 --sample-expressions Sap68 --sample-expressions SRR4162919 --sample-expressions SRR4162934 --sample-expressions SRR4235992 --sample-expressions SRR4236033 --sample-expressions Sul2 --sample-expressions Sul4 --sample-expressions Sul5 --sample-expressions Tar12 --sample-expressions Tar1 --exclude-sample-expressions D50 --exclude-sample-expressions D51 --exclude-sample-expressions D52 --exclude-sample-expressions D53 --exclude-sample-expressions D54 --exclude-sample-expressions D55 --exclude-sample-expressions D56 --exclude-sample-expressions D57 --exclude-sample-expressions D58 --exclude-sample-expressions D59 --exclude-non-variants true -O ${NAME}_A_Italia.vcf
gatk SelectVariants -R ./Bdistachyon_556_v3.0.fa -V $i --exclude-sample-expressions D71 --exclude-sample-expressions D72 --exclude-sample-expressions Cb23 --exclude-sample-expressions Cro24 --exclude-sample-expressions D65 --exclude-sample-expressions D66 --exclude-sample-expressions D67 --exclude-sample-expressions D68 --exclude-sample-expressions D69 --exclude-sample-expressions D70 --exclude-sample-expressions Mca12 --exclude-sample-expressions Mca16 --exclude-sample-expressions Mca18 --exclude-sample-expressions Mca19 --exclude-sample-expressions Mca23 --exclude-sample-expressions D25 --exclude-sample-expressions Msa11 --exclude-sample-expressions Msa16 --exclude-sample-expressions Msa6 --exclude-sample-expressions Msa9 --exclude-sample-expressions San11 --exclude-sample-expressions San12 --exclude-sample-expressions San25 --exclude-sample-expressions Spm16 --exclude-sample-expressions SRR4162922 --exclude-sample-expressions SRR4162932 --exclude-sample-expressions SRR4162933 --exclude-non-variants true -O ${NAME}_full_A_B_cluster.vcf
done

## the chromosomes in the vcf need to be numberd as follow: 1,2,3,4,5
for i in  $(ls Phased_Pol_SNPs_TE_Bd*_*.vcf)
do 
sed -i -e 's/Bd1/1/g' $i
sed -i -e 's/Bd2/2/g' $i
sed -i -e 's/Bd3/3/g' $i
sed -i -e 's/Bd4/4/g' $i
sed -i -e 's/Bd5/5/g' $i
done

for i in  $(ls Phased_Pol_SNPs_TE_Bd*_*.vcf)
do 
gzip $i
done

## run geva:
## example using the A_East clade
## data convertion GEVA:

./GEVA/geva-master/geva_v1beta --vcf Phased_Pol_SNPs_TE_Bd1_A_East.vcf.gz --map ./Bdistachyon_genetic_map_chr1.txt --out Bd1_A_East
./GEVA/geva-master/geva_v1beta --vcf Phased_Pol_SNPs_TE_Bd2_A_East.vcf.gz --map ./Bdistachyon_genetic_map_chr2.txt --out Bd2_A_East
./GEVA/geva-master/geva_v1beta --vcf Phased_Pol_SNPs_TE_Bd3_A_East.vcf.gz --map ./Bdistachyon_genetic_map_chr3.txt --out Bd3_A_East
./GEVA/geva-master/geva_v1beta --vcf Phased_Pol_SNPs_TE_Bd4_A_East.vcf.gz --map ./Bdistachyon_genetic_map_chr4.txt --out Bd4_A_East
./GEVA/geva-master/geva_v1beta --vcf Phased_Pol_SNPs_TE_Bd5_A_East.vcf.gz --map ./Bdistachyon_genetic_map_chr5.txt --out Bd5_A_East

## make a list of the site possitions and remove first and last position as they cause problems
zcat Phased_Pol_SNPs_TE_Bd1_A_East.vcf.gz | cut -f 2 | sed '1,62d' | sed '$ d' > Bd1_A_East_pos_list
zcat Phased_Pol_SNPs_TE_Bd2_A_East.vcf.gz | cut -f 2 | sed '1,62d' | sed '$ d' > Bd2_A_East_pos_list
zcat Phased_Pol_SNPs_TE_Bd3_A_East.vcf.gz | cut -f 2 | sed '1,62d' | sed '$ d' > Bd3_A_East_pos_list
zcat Phased_Pol_SNPs_TE_Bd4_A_East.vcf.gz | cut -f 2 | sed '1,62d' | sed '$ d' > Bd4_A_East_pos_list
zcat Phased_Pol_SNPs_TE_Bd5_A_East.vcf.gz | cut -f 2 | sed '1,62d' | sed '$ d' > Bd5_A_East_pos_list

## Split the position files in smaller files to be able to run GEVA more efficiently:
for i in {1..167}
do
START=`expr  \( $i - 1 \) \* 5000 + 1`
END=`expr 5000 \* $i`
sed -n "${START},${END}p"  Bd1_A_East_pos_list > Bd1_A_East_pos_list_${START}_${END}
done

for i in {1..143}
do
START=`expr  \( $i - 1 \) \* 5000 + 1`
END=`expr 5000 \* $i`
sed -n "${START},${END}p"  Bd2_A_East_pos_list > Bd2_A_East_pos_list_${START}_${END}
done

for i in {1..153}
do
START=`expr  \( $i - 1 \) \* 5000 + 1`
END=`expr 5000 \* $i`
sed -n "${START},${END}p"  Bd3_A_East_pos_list > Bd3_A_East_pos_list_${START}_${END}
done

for i in {1..106}
do
START=`expr  \( $i - 1 \) \* 5000 + 1`
END=`expr 5000 \* $i`
sed -n "${START},${END}p"  Bd4_A_East_pos_list > Bd4_A_East_pos_list_${START}_${END}
done

for i in {1..69}
do
START=`expr  \( $i - 1 \) \* 5000 + 1`
END=`expr 5000 \* $i`
sed -n "${START},${END}p"  Bd5_A_East_pos_list > Bd5_A_East_pos_list_${START}_${END}
done

## estimate age
## using the Run_GEVA_age_estimate.lsf script
ls Bd1_A_East_pos_list_* > GEVA_list_of_position_files.txt
ls Bd2_A_East_pos_list_* >> GEVA_list_of_position_files.txt
ls Bd3_A_East_pos_list_* >> GEVA_list_of_position_files.txt
ls Bd4_A_East_pos_list_* >> GEVA_list_of_position_files.txt
ls Bd5_A_East_pos_list_* >> GEVA_list_of_position_files.txt

bsub < ./Run_GEVA_age_estimate.lsf -J "GEVA[1-638]"

## merge result files:
head -1 Bd1_A_East_400001_405000.sites2.txt > header.sites2.txt

for j in  Bd1 Bd2 Bd3 Bd4 Bd5
do

for i in  $(ls header.sites2.txt ${j}_A_East_*0.sites2.txt | sort -r)
do
if [ $i == "header.sites2.txt" ]
then
cat $i
else
grep " J " $i
fi
done > ${j}_A_East_Clock_J_age_estimate.sites2.txt

for i in  $(ls header.sites2.txt ${j}_A_East_*0.sites2.txt | sort -r)
do
if [ $i == "header.sites2.txt" ]
then
cat $i
else
grep " M " $i
fi
done > ${j}_A_East_Clock_M_age_estimate.sites2.txt

for i in  $(ls header.sites2.txt ${j}_A_East_*0.sites2.txt | sort -r)
do
if [ $i == "header.sites2.txt" ]
then
cat $i
else
grep " R " $i
fi
done > ${j}_A_East_Clock_R_age_estimate.sites2.txt

done




