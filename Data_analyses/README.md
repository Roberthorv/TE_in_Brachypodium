The script Age-adjusted_SFS.R can be used to generate the age adjusted SFS for retrotransposons.
The inputs needed are the GEVA age estimates for TEs and SNPs, a vcf with the synonymous sites, a table with the IDs of retrotransposons (can be generated using the Make_TIP_TAP_info_file.R script). 

The age adjusted SFS for retrotransposons split by distence to the next gene can be done using the Age-adjusted_SFS_based_on_gene_distance.R
The age adjusted SFS for DNA-transposons split by distence to the next gene can be done using the Age-adjusted_SFS_DNA_transposons.R

The PCA, ANCOVA and iHS analyses can be done using the Run_PCA.R, Run_ANCOVA.R and iHS_analyses_combine_regions.sh script 
