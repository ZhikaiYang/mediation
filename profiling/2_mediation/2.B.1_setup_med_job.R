### Jinliang Yang
### 01-22-2019
### set up array job

##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

JOBID <- as.numeric(as.character(args[1]))
print(JOBID)

library(glmnet)
library(MASS)
library(rrBLUP)
library(parallel)
library(doParallel)
library(CompQuadForm)

library(dplyr)
library("data.table")

source('lib/highmed2019.r')
source('lib/fromSKAT.R')
source('lib/MedWrapper.R')
source('lib/reporters.R')
source('lib/slurm_wrapper.R')

df <- read.csv("largedata/df_job_control.csv")

slurm_med_wrapper(ti=as.numeric(as.character(df$trait[JOBID])), cutoff_pm=0.05, ncores=8,
                  phenofile = "data/geno_trait.txt",
                  genofile = "largedata/geno/gddtosilk_p1msnps_matrix.txt",
                  rnafile = as.character(df$rnafile[JOBID]),
                  pcfile="largedata/allchr_bisnp_n282_snpid_maf01_geno2_pruned.eigenvec",
                  outdir="largedata/med_output")

#####

