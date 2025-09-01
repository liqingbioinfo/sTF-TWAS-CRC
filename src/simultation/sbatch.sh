#!/bin/bash
#SBATCH --job-name=sim_y
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0:0:0

#conda activate tftwas_coloc
Rscript 1.simulate_Yvalues_basedonGeno.R $1
#Rscript 3.glm_association.R $1
