#!/bin/bash
#SBATCH --job-name=glm_logi
#SBATCH --output=glm_logi%A_%a.out
#SBATCH --error=glm_logi%A_%a.err
#SBATCH --time=7-0:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8

# 1000G_EUR_Phase3_plink_hg38 download from https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays 
geno_file="/ref/Alkesgroup/1000G_EUR_Phase3_plink_hg38/merged/1000G.EUR.QC.1_22-merge"
pheno_file="/CRC_TF_TWAS_zhishan2024/NC_revision_OceanCode/v2_NC_revision/simulation_values_res/sim_y_values/beta1_1_h2_0.5_SNPs_1000.sim_y"

read -a list < <(awk '{for(i=9;i<=NF;++i) printf "%s ", $i; print ""}' ${pheno_file}.psam | head -n1)
for item in "${list[@]}"; do
    echo "$item"
    plink2 --pfile ${geno_file} --pheno ${pheno_file}.psam --pheno-name $item --covar ${pheno_file}.psam --covar-name SEX,PC1,PC2,PC3,PC4,PC5 --glm hide-covar --out ${pheno_file}/10M
done
