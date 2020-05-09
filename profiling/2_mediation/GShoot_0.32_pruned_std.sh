#!/bin/sh
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=150gb
#SBATCH --time=12:00:00
#SBATCH --partition=batch
#SBATCH --licenses=common
#SBATCH --job-name=GShoot_0.032_pruned_std
#SBATCH --mail-user=zhikaiyang911@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/GShoot_0.032_pruned_std.err
#SBATCH --output=/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/GShoot_0.032_pruned_std.out
module load R/3.5
Rscript /common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/GShoot_0.032_pruned_std.R
