#!/bin/bash
#SBATCH --job-name=ldscore_chr
#SBATCH --output=ldscore_chr%A_%a.out
#SBATCH --error=ldscore_chr%A_%a.err
#SBATCH --time=7-0:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

# 1000G_EUR_Phase3_plink_hg38 download from https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays 
PREF="Your_folder_name"
PLINK_PREFIX="$PREF/ref/Alkesgroup/1000G_EUR_Phase3_plink_hg38"
WK_DIR="${PREF}/84Tracks/"
LIST="${WK_DIR}/84.tracks.list.txt"

###Make annotation
while read -r FILE; do
	#make annot file
	for chr in {1..22}; do
		${PREF}/biosoft/ldsc/make_annot.py \
		--bed-file "${FILE}_peaks.bed.bed" \
		--bimfile "${PLINK_PREFIX}/1000G.EUR.QC.${chr}.bim" \
		--annot-file "${WK_DIR}/annots/${FILE}.${chr}.annot.txt"
	done
done < "$LIST"

###Merge multiple annotation into one for each chrom

###Make ldscores
ANNOT_DIR="${PREF}/84Tracks/merge_tracks_annots/"
OUTPUT_DIR="${PREF}/84Tracks/ldscores/"

mkdir -p "$OUTPUT_DIR"

for CHR in {1..22}; do
	python ${PREF}/biosoft/ldsc/ldsc.py \
		--l2 \
		--bfile ${PREF}/ref/Alkesgroup/1000G_EUR_Phase3_plink_hg38/1000G.EUR.QC.${CHR} \
		--ld-wind-kb 100 \
		--annot ${ANNOT_DIR}/${CHR}.annot.txt \
		--out ${OUTPUT_DIR}/${CHR}
done