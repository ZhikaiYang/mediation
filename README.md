# The code for this project is in profiling folder;

inside of it, there are 1_GWAS folder, which holds the code for doing GWAS:

The prefix of the code name indicates the order to be carried out;
Codes needs to be done in bash environment is in bash chunks;
Codes needs to be done in R environment is in R chunks, if it is needed to be done in high performance computing environment, it will be noted in the front of the R chunks.
To carry out the analysis, some part of the code needs to be adjusted like the path.
For simplicity, the result files needed for manhattan plot for GWAS are moved into data folder.


# Results

### Meidation analysis results

- direct SNP `data/input/dsnps_20149rows.csv`
- indirect SNP `data/input/isnps_55257rows.csv`
- mediator genes `data/input/mediators_23296rows_by_cat.csv`

### GWAS results
raw data:
- agronomic triats GWAS :`/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/geno_282/output`
- metabolomics triats GWAS : `/common/jyanglab/shared/Gen_Xu/282_metabolite/output`

processed data:

-all traits GWAS in merged significant SNPs region within 100kb: `/common/jyanglab/zhikaiyang/projects/mediation/largedata/dsnps_vs_gwas/all_traits.all_res.txt`

# For Mediation analysis, the code is in 2_mediation folder;

there is no specific order for each R code, and they were run in high performance computing environment, and the jobs are submited with their corresponding slurm files;
To carry out the analysis, some part of the code needs to be adjusted like the path.
For simplicity, the result files needed for ploting the mediator genes in manhattan plot are moved into data folder.


For visualization the mediation analysis result in manhattan plot, go the 3_visualization folder;
The part of code that is needed to be done in HCC are noted;
And for simplicity, the result generated in HCC are moved into data folder

The plots in Fig.6 in the paper are from the two commands with comment "#mediator genes ploted in manhattan plot" above for the two tissues germinating shoot and leaf tips, noted by the comment; these two plots are each the first plot for each tissues in the code
The 2nd plot for each tissue is the plot showing the position of each direct snps(snps whose coefficients are significantly non-zero in outcome model)
The 3rd plot for each tissue is the plot showing the position of each indirect snps(snps whose coefficients are significantly non-zero in mediator model) 




