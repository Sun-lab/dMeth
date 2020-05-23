#!/bin/bash
#SBATCH --job-name runLUSC 
#SBATCH --mail-user=zhanyu@fredhutch.org
#SBATCH --mail-type=ALL
#SBATCH --nodes=1 #SBATCH --output=Rout/par-%J.out
#SBATCH --error=Rout/par-%J.err
#SBATCH --cpus-per-task=1
ml R/3.6.0-foss-2016b-fh1
R CMD BATCH runLUSC.R Rout/runLUSC.Rout
