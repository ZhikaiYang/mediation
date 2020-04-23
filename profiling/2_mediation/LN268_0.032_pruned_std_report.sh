#!/bin/sh
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=60gb
#SBATCH --time=12:00:00
#SBATCH --partition=jyanglab
#SBATCH --licenses=common
#SBATCH --job-name=LN268_0.032_pruned_std_report
#SBATCH --mail-user=zhikaiyang911@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/LN268_0.032_pruned_std_report.err
#SBATCH --output=/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/LN268_0.032_pruned_std_report.out
module load R/3.5
Rscript /common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/LN268_0.032_pruned_std_report.R
