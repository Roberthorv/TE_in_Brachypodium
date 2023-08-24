# example of a loop to run SLiM with: negative selection coefficient of -0.025, proportion of netral mutation 0% and 80% selfing

for ((i=1; i<=20; i++))
do
echo $i
./SLiM/build/slim -d run=$i -d "name='SNPs_sel_-0.025_proportion_1'" -d negative=-0.025 -d proportion=1 -d selfing=0.80 B.distachyon_pop_scaling_50_SNPs_for_shushu_selfing_var.slim
done

