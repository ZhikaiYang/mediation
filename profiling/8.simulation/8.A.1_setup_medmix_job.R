### Jinliang Yang
### 12-22-2021
### set up array job

##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

JOBID <- as.numeric(as.character(args[1]))
print(JOBID)

library(data.table)
library(glmnet)
library(MASS)
library(rrBLUP)
library(parallel)
library(doParallel)
library(CompQuadForm)
source('lib/highmed2019.r')
source('lib/fromSKAT.R')
source("profiling/8.simulation/8.A.1_setup_medmix_job.R")

df <- read.csv("profiling/8.simulation/df_job_control_medmix.csv")

mediation_run(nQTLperM=as.numeric(as.character(df$nQTLper[JOBID])), 
              nQTLar=as.numeric(as.character(df$nQTLar[JOBID])), 
              h2_as=as.numeric(as.character(df$h2_as[JOBID])), 
              seedas.numeric(as.character(df$seed[JOBID])), 
              ncores = 4)
#####

