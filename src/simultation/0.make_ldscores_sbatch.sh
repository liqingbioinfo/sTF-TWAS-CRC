#!/bin/bash
#SBATCH --job-name=ldscore_chr
#SBATCH --output=ldscore_chr%A_%a.out
#SBATCH --error=ldscore_chr%A_%a.err
#SBATCH --time=7-0:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

# 1000G_EUR_Phase3_plink_hg38 download from https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays 
PREF="/data/sbcs/GuoLab/backup/liq17/CRC_TF_TWAS_zhishan2024/NC_revision_OceanCode/NC_revision"
PLINK_PREFIX="/data/sbcs/GuoLab/backup/liq17/ref/Alkesgroup/1000G_EUR_Phase3_plink_hg38"
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

###Make ldscores [multi annot per chr, 10M]
ANNOT_DIR="${PREF}/84Tracks/multi_annots_ldscores"
OUTPUT_DIR="${PREF}/84Tracks/multi_annots_ldscores/ldscores_10M/"
mkdir -p "$OUTPUT_DIR"
for CHR in {1..22}; do
	python /data/sbcs/GuoLab/backup/liq17/biosoft/ldsc/ldsc.py \
		--l2 \
		--bfile ${PLINK_PREFIX}/1000G.EUR.QC.${CHR} \
		--ld-wind-kb 100 \
		--annot ${ANNOT_DIR}/${FILE}.${CHR}.annot.txt \
		--out ${OUTPUT_DIR}/${CHR}
done

###Make ldscores [per annot per chr], Not recommanded, slow
while read -r FILE; do
	ANNOT_DIR="${PREF}/84Tracks/one_annot_ldscores/${FILE}"
	OUTPUT_DIR="${PREF}/84Tracks/one_annot_ldscores/ldscores/${FILE}"
	mkdir -p "$OUTPUT_DIR"
	for CHR in {1..22}; do
		python /data/sbcs/GuoLab/backup/liq17/biosoft/ldsc/ldsc.py \
			--l2 \
			--bfile ${PLINK_PREFIX}/1000G.EUR.QC.${CHR} \
			--ld-wind-kb 100 \
			--annot ${ANNOT_DIR}/${FILE}.${CHR}.annot.txt --thin-annot  \
			--out ${OUTPUT_DIR}/${CHR}
	done
done < "$LIST"


###Make ldscores [multi annot per chr, 100K]

###Make annotation
while read -r FILE; do
	#make annot file
	for chr in {1..22}; do
		/data/sbcs/GuoLab/backup/liq17/biosoft/ldsc/make_annot.py \
		--bed-file "${FILE}_peaks.bed.bed" \
		--bimfile "${PLINK_PREFIX}/processed/${chr}.nodup.sorted.maf01.LD50_5_01.1000G_hg38.filtered.bim" \
		--annot-file "${WK_DIR}/multi_annots_ldscores/annot_100K/${FILE}.${chr}.annot.txt"
	done
done < "$LIST"

ANNOT_DIR="${PREF}/84Tracks/multi_annots_ldscores/annot_100K/"
OUTPUT_DIR="${PREF}/84Tracks/multi_annots_ldscores/ldscores_100K/"
mkdir -p "$OUTPUT_DIR"
for CHR in {1..22}; do
	python /data/sbcs/GuoLab/backup/liq17/biosoft/ldsc/ldsc.py \
		--l2 \
		--bfile ${PLINK_PREFIX}/processed/${CHR}.nodup.sorted.maf01.LD50_5_01.1000G_hg38.filtered \
		--ld-wind-kb 100 \
		--annot ${ANNOT_DIR}/${CHR}.annot.txt \
		--out ${OUTPUT_DIR}/${CHR}
done


