## Indexing ref:
# bowtie2:
bowtie2-build Bdistachyon_556_v3.0.fa Bdistachyon_556_v3.0

# yaha:
/media/HDD2/Robert/yaha/bin/yaha -g Bdistachyon_556_v3.0.fa

## Formating the TE annotation file:
wc -l TE_annotation/Bdistachyon_314_v3.0.fa.mod.EDTA.TEanno.gff3
tail -196274 TE_annotation/Bdistachyon_314_v3.0.fa.mod.EDTA.TEanno.gff3 > int_a
cut -f 1,4,5,7 int_a > part_a
cut -f 9 int_a > int_b
sed -e 's/;Name=/\t/g' int_b | cut -f 2 | sed -e 's/;/\t/g' | cut -f 1 > part_b
cut -f 8 int_a > part_c
sed -e 's/;Classification=/\t/g' int_b | cut -f 2 | sed -e 's/;/\t/g' | cut -f 1 > part_d

paste part_a part_b part_c part_d > Bdistachyon_314_v3.0.fa.mod.EDTA.TEanno.bed

rm int_a int_b part_a part_b part_c part_d
