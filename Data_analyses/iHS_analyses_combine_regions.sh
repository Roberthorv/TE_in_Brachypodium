## combine all regions for the ANCOVA

sed "s/$/\tB_West/" high_iHS_regions_B_West_pos_sel_less_5_perc.bed > high_iHS_regions_B_West_pos_sel_less_5_perc_ID_tab.bed
sed "s/$/\tB_East/" high_iHS_regions_B_East_pos_sel_less_5_perc.bed > high_iHS_regions_B_East_pos_sel_less_5_perc_ID_tab.bed
sed "s/$/\tA_East/" high_iHS_regions_A_East_pos_sel_less_5_perc.bed > high_iHS_regions_A_East_pos_sel_less_5_perc_ID_tab.bed
sed "s/$/\tA_Italy/" high_iHS_regions_A_Italy_pos_sel_less_5_perc.bed > high_iHS_regions_A_Italy_pos_sel_less_5_perc_ID_tab.bed
cat high_iHS_regions_B_West_pos_sel_less_5_perc_ID_tab.bed high_iHS_regions_B_East_pos_sel_less_5_perc_ID_tab.bed high_iHS_regions_A_East_pos_sel_less_5_perc_ID_tab.bed high_iHS_regions_A_Italy_pos_sel_less_5_perc_ID_tab.bed | bedtools sort | bedtools merge -i - -c 4 -o collapse > high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab.bed
bedtools intersect -a high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab.bed -b Full_clean_TIP_TAP_Sample_table_326_samples_simlpe.bed -wb > TEs_in_high_iHS_regions_all_clades_pos_sel_less_5_perc.bed
