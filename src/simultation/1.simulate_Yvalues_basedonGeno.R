#########################################################
####  Computer system requirements                    ###
####  R package:  data.table, lme4                    ###
####  Unix operating system environment               ###
####  20 GB memory for per_geno, 120G for whole geno  ### 
#########################################################
# conda activate tftwas
library(data.table) # read in large data matrix fast
library(lme4)       # for mixed model
library(parallel)
library(dplyr)
args = commandArgs(trailingOnly = TRUE)
setwd("/CRC_TF_TWAS_zhishan2024/NC_revision_OceanCode/v2_NC_revision/")

### Load TF annotation matrix
TFs_annots<-as.data.frame(fread(paste0("./84Tracks/annots/1_22.annot.txt")))
# remove duplicate variants 
TFs_annots <- TFs_annots[!duplicated(TFs_annots$SNP),]
TFs_annots<-TFs_annots[TFs_annots$CHR>=1 & TFs_annots$CHR<=22,]

### Select causal SNPs and make causal SNP annotation files
set.seed(20250828) 
n_causal=as.numeric(args[1])
n_round=50
beta_1=1
h2=0.5
SNP_causal_pool = TFs_annots[TFs_annots$ANNO==1,]$SNP
sim_y_matrix = matrix(integer(0), nrow = n_round , ncol= 489) %>% as.data.frame()

# start from TFs_annots with CHR,BP,SNP,CM
DT <- as.data.table(TFs_annots[, c("CHR","BP","SNP","CM")])
simNames <- paste0("SIM", seq_len(n_round))
# store each round’s sampled SNPs for validation
causal_list <- vector("list", n_round)
for (sim_i in seq_len(n_round)) {
  SNP_causal <- sample(SNP_causal_pool, n_causal)
  causal_list[[sim_i]] <- SNP_causal
  nm <- simNames[sim_i]
  DT[, (nm) := 0L]
  DT[data.table(SNP = SNP_causal), on = .(SNP), (nm) := 1L]
}
# back to data.frame (optional)
SNP_causal_annots <- as.data.frame(DT)
#Validate SNPs: all( (SNP_causal_annots$SNP %in% causal_list[[2]]) == (SNP_causal_annots$SIM2 == 1) )

for(chr in 1:22){
	SNP_causal_annots_per_chr = SNP_causal_annots[SNP_causal_annots$CHR==as.character(chr), ]
	print(dim(SNP_causal_annots_per_chr))
	fwrite(SNP_causal_annots_per_chr, paste0("./84Tracks/annots_sim", n_causal,"/",chr,".annot.txt") ,sep="\t")
}

for(sim_i in 1:n_round){
	simName=paste0("SIM", sim_i)
	# >> make ldscores absed on the simulated annotations
	### Simulate Y based on causal SNPs
	SNP_causal_geno=data.frame()
	for(chr in 1:22){
		SNP_causal_chr = SNP_causal_annots[(SNP_causal_annots[simName]==1) & (SNP_causal_annots$CHR == as.character(chr)), ]$SNP
		print(paste0("Num of annotated SNPs is for ", chr, " is ",length(SNP_causal_chr)))
		chr_geno=data.frame(fread(paste0("/ref/Alkesgroup/1000G_EUR_Phase3_plink_hg38/processed/",chr,".nodup.sorted.1000G_hg38.traw")))
		SNP_causal_chr_geno = chr_geno[chr_geno$SNP %in% SNP_causal_chr, ]
		rownames(SNP_causal_chr_geno) <- SNP_causal_chr_geno$SNP
		SNP_causal_chr_geno_clean <- SNP_causal_chr_geno[, 7:dim(SNP_causal_chr_geno)[2]] #select only genotype
		SNP_causal_geno=rbind(SNP_causal_geno, SNP_causal_chr_geno_clean)
		rm(chr_geno)
	}
	print(paste0("Total num of annotated SNPs is ",dim(SNP_causal_geno)[1]))
	SNP_causal_geno_r0a1 <- 2 - SNP_causal_geno  ##Chek the coding for ref and alt from *.traw and *.tped file
	tmp = sub("_(.*)", "", colnames(SNP_causal_geno_r0a1))
	colnames(SNP_causal_geno_r0a1) <- tmp
	sample_size= dim(SNP_causal_geno_r0a1)[2]
	y_geno <-colSums(beta_1 * SNP_causal_geno_r0a1)
	y_err <- rnorm(sample_size, mean=0, sd=sqrt((1-h2)/h2*var(y_geno)))
	y_sim = y_geno + y_err
	#y_sim=y_geno
	y_sim_cc <- ifelse(y_sim >= median(y_sim), 1, 2) #1 for case, 2 for control
	sim_y_matrix[sim_i, ]= y_sim_cc
}
sim_y_matrix_t = transpose(sim_y_matrix)
sim_y_df_t <- as.data.frame(sim_y_matrix_t)  #sample N x simulate times t
sim_y_df_t["IID"]=names(y_sim_cc)  # Add IID column
sim_y_df_t_IID = setcolorder(sim_y_df_t, c("IID", setdiff(names(sim_y_df_t), "IID"))) # reorder col
colnames(sim_y_df_t_IID) <- c("IID", simNames)
fwrite(sim_y_df_t_IID,paste0("./simulation_values_res/sim_y_values/beta1_1_h2_",h2,"_SNPs_",n_causal,".sim_y.txt"), sep="\t", row.names = FALSE, col.names = TRUE)
