## Prepare the genetic map called Bdistachyon_genetic_map.txt
sed -i -e 's/_//g' Bdistachyon_genetic_map.txt

## split genetic map by scaffolds:
grep "chr1\|Chromosome" Bdistachyon_genetic_map.txt > Bdistachyon_genetic_map_chr1.txt
grep "chr2\|Chromosome" Bdistachyon_genetic_map.txt > Bdistachyon_genetic_map_chr2.txt
grep "chr3\|Chromosome" Bdistachyon_genetic_map.txt > Bdistachyon_genetic_map_chr3.txt
grep "chr4\|Chromosome" Bdistachyon_genetic_map.txt > Bdistachyon_genetic_map_chr4.txt
grep "chr5\|Chromosome" Bdistachyon_genetic_map.txt > Bdistachyon_genetic_map_chr5.txt


