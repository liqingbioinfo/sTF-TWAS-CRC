#!/bin/bash
#SBATCH --job-name=ldscore_chr
#SBATCH --output=ldscore_chr%A_%a.out
#SBATCH --error=ldscore_chr%A_%a.err
#SBATCH --time=7-0:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

# 1000G_EUR_Phase3_plink_hg38 download from https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays 
PREF="/CRC_TF_TWAS_zhishan2024/NC_revision_OceanCode/v2_NC_revision/"
PLINK_PREFIX="/ref/Alkesgroup/1000G_EUR_Phase3_plink_hg38"
WK_DIR="${PREF}/84Tracks/"


for chr in {1..22}; do
	/biosoft/ldsc/make_annot.py \
	--bed-file "${WK_DIR}/merged_all.sorted.bed" \
	--bimfile "${PLINK_PREFIX}/1000G.EUR.QC.${chr}.bim" \
	--annot-file "${WK_DIR}/annots/${chr}.annot.txt"
done

##ldscores based on annotations from 600K occupied SNPs from 84 tracks
ANNOT_DIR="${WK_DIR}/annots/"
OUTPUT_DIR="${WK_DIR}/ldscores/"
mkdir -p "$OUTPUT_DIR"
for CHR in {1..22}; do
	python /biosoft/ldsc/ldsc.py \
		--l2 \
		--bfile ${PLINK_PREFIX}/1000G.EUR.QC.${CHR} \
		--ld-wind-kb 100 \
		--annot ${ANNOT_DIR}/${CHR}.annot.txt --thin-annot\
		--out ${OUTPUT_DIR}/${CHR}
done

##ldscores based on annotations from simulated annotations

ANNOT_DIR="${WK_DIR}/annots_sim1000/"
OUTPUT_DIR="${WK_DIR}/ldscores_sim1000/"
mkdir -p "$OUTPUT_DIR"
for CHR in {11..22}; do
	python /biosoft/ldsc/ldsc.py \
		--l2 \
		--bfile ${PLINK_PREFIX}/1000G.EUR.QC.${CHR} \
		--ld-wind-kb 100 \
		--annot ${ANNOT_DIR}/${CHR}.annot.txt \
		--out ${OUTPUT_DIR}/${CHR}
done





