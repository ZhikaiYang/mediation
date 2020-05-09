#!/bin/sh
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=150gb
#SBATCH --time=12:00:00
#SBATCH --partition=batch
#SBATCH --licenses=common
#SBATCH --job-name=kern_0.032_pruned_std_report
#SBATCH --mail-user=zhikaiyang911@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/kern_0.032_pruned_std_report.err
#SBATCH --output=/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/kern_0.032_pruned_std_report.out
module load R/3.5
Rscript /common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/kern_0.032_pruned_std_report.R
