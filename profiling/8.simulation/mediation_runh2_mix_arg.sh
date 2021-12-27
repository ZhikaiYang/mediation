#!/bin/sh
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --time=90:00:00
#SBATCH --partition=batch
#SBATCH --licenses=common
#SBATCH --job-name=mediation_runh2_mix
#SBATCH --mail-user=zhikaiyang911@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -o /common/jyanglab/zhikaiyang/projects/mediation/slurm-log/mediation_runh2_mix.txt
#SBATCH -e /common/jyanglab/zhikaiyang/projects/mediation/slurm-log/mediation_runh2_mix.txt


module load R/3.5
Rscript /common/jyanglab/zhikaiyang/projects/mediation/largedata/simulation/simulation_corn_data_mediation_test_h2_mix.R  nQTLperM nQTLar h2_as seed

