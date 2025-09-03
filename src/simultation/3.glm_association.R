##################
# conda activate tftwas_coloc
# 80G
##################

library(data.table)
library(lme4)
library(dplyr)
args = commandArgs(trailingOnly = TRUE)

set.seed(123)
KB<-100000
glm_res<-matrix(integer(0), nrow=100, ncol=6) %>% as.data.frame()
names(glm_res) <- c("Beta1", "CistromeID", "Slope", "SD", "Tvalue", "Pvalue")
wk="/data/sbcs/GuoLab/backup/liq17/CRC_TF_TWAS_zhishan2024/NC_revision_OceanCode/v2_NC_revision"
#beta1=args[1]
beta1=1
h2=0.5
snps=as.character(args[1])

###load sim_ys
y_sim=data.frame(fread(paste0(wk, "/simulation_values_res/sim_y_values/beta1_",beta1,"_h2_",h2,"_SNPs_",snps,".sim_y.txt"), header = TRUE))
tmp=colnames(y_sim)[2:length(colnames(y_sim))]
tracks_num = sub("^X", "", tmp)

###load annots
annots=data.frame()
for(chr in 1:22){
	annot_per_chr=fread(paste0(wk,"/84Tracks/annots_sim",snps,"/",chr,".annot.txt"))
	annots=rbind(annots, annot_per_chr)
}
annots=as.data.frame(annots)

#y_track = colnames(y_sim)[3]
#y_track_name=paste0("X", y_track)
count=1
for(y_track in tracks_num){
	y_track_name=y_track
	gwas_file=paste0(wk,"/simulation_values_res/sim_y_values/beta1_",beta1,"_h2_",h2,"_SNPs_",snps,".sim_y/10M.",y_track,".glm.logistic.hybrid")
	###load GWAS SS
	if(file.exists(gwas_file)){
		gwas_ss=fread(gwas_file)
		names(gwas_ss)[names(gwas_ss)=="ID"]<-"SNP"
		gwas_ss <- gwas_ss[!is.na(Z_STAT) & `#CHROM` >= 1 & `#CHROM` <= 22]  #Z_STAT  P
		gwas_ss$A2 = ifelse(gwas_ss$A1==gwas_ss$ALT, gwas_ss$REF, gwas_ss$ALT)
		gwas_ss$N=489
		ldsc_gwas_ss=gwas_ss[,c("SNP", "N", "A1", "A2", "Z_STAT", "P")]
		colnames(ldsc_gwas_ss) <- c("SNP", "N", "A1", "A2", "Z", "P")
		fwrite(ldsc_gwas_ss, paste0(wk,"/simulation_values_res/sim_y_values/beta1_",beta1,"_h2_",h2,"_SNPs_",snps,".sim_y/10M.",y_track,".glm.logistic.ldsc.txt"), sep="\t")
		
		###merge GWAS SS and Annots
		gwas_ss_annots=left_join(gwas_ss, annots[, c("SNP", y_track_name)], by = "SNP")
		gwas_ss_annots$loci <- paste0(gwas_ss_annots[["#CHROM"]],'_',floor(gwas_ss_annots[["POS"]]/KB)) # CHR and POS are from annots files

		###Run glm and output
		out <- summary(lmer(I(gwas_ss_annots$Z_STAT^2) ~ gwas_ss_annots[[y_track_name]]+(1|gwas_ss_annots$loci),control = lmerControl(calc.derivs = FALSE)))
		####One-tail p value 
		if(out$coef[2,1]<0){
			p_ <- pnorm(out$coef[2,3])
		}else{
			p_ <- 1 - pnorm(out$coef[2,3])
		}
		print(c(beta1, y_track, out$coef[2,1], out$coef[2,2], out$coef[2,3], p_))
		glm_res[count,]=c(beta1, y_track, out$coef[2,1], out$coef[2,2], out$coef[2,3], p_)
		count=count+1
	}
}
fwrite(as.data.frame(glm_res), paste0(wk, "/simulation_values_res/lmer_res/beta1_",beta1,"_h2_",h2,"_SNPs_",snps,".sim_y.lmer.csv"))
