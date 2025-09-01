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
setwd("/data/sbcs/GuoLab/backup/liq17/CRC_TF_TWAS_zhishan2024/NC_revision_OceanCode/NC_revision")

# Input TF annotation matrix
TFs_annot <- as.data.frame(fread("./data/1_22.annot.txt",header=T))
names(TFs_annot)
## remove duplicate variants 
TFs_annot <- TFs_annot[!duplicated(TFs_annot$SNP),] #MarkerName	Allele1	Allele2	Freq1	FreqSE	MinFreq	MaxFreq	Effect	StdErr	P-value	Direction	HetISq	HetChiSq	HetDf	HetPVal	TotalSampleSize
TFs_annot<-TFs_annot[TFs_annot$CHR>=1 & TFs_annot$CHR<=22,]

# load Genotype [RAM 120G OK, 80G killed]
#geno_raw <- fread("./data/1000G_ref_hg38.traw") # Ref coded as 2 and alt coded as 0
#geno_raw_with_annot <- geno_raw[rowSums(geno_raw[ , 7:ncol(geno_raw)]) != 0, ] # Select all rows with at least 1 annot
#rownames(geno_raw_with_annot) <- geno_raw_with_annot$SNP
#geno_clean =2 - geno_raw_with_annot[ , 7:ncol(geno_raw_with_annot)] # swith ref and alt code
#rm(geno_raw)
#rm(geno_raw_with_annot)
#tmp = sub("_(.*)", "", colnames(geno_clean))
#colnames(geno_clean) <- tmp

#### simulate Y_ij based on beta_1 and beta_0
#### Y_ij = beta_1 x X_ij + episilon, beta_1 from a list and episilong ~ N(0, 1)
tf_cols=84
n_round=60
h2s=c(0.05, 0.3, 0.6)
h2 = h2s[as.numeric(args[1])]
beta_1 = 1
set.seed(250822)
sim_y_list = list()
sim_track_name = vector(mode="character", n_round)
lmer_res_list = list()
sim_y_matrix = matrix(integer(0), nrow = n_round , ncol= 489) %>% as.data.frame()

for (sim_i in 1:n_round) {
  print(paste0("sim round " , as.character(sim_i), "/", as.character(n_round)))
  #sim_i_track_index = sample(5:dim(TFs_annot)[2], 1)
  #sim_i_track_name = colnames(TFs_annot)[sim_i_track_index]
  #sim_track_name[sim_i] <-sim_i_track_name #record track name
  sim_i_track_name <- colnames(TFs_annot)[sim_i+4]
  sim_track_name[sim_i] <- sim_i_track_name
  geno_annot_clean=data.frame()
  for(chr in 1:22){
    chr_annot = TFs_annot[(TFs_annot$CHR==chr) & (TFs_annot[[sim_i_track_name]]==1), ]$SNP
    print(paste0(sim_i_track_name, " : Num of annotated SNPs is ",length(chr_annot), " for ", chr))
    chr_geno=data.frame(fread(paste0("/data/sbcs/GuoLab/backup/liq17/ref/Alkesgroup/1000G_EUR_Phase3_plink_hg38/processed/",chr,".nodup.sorted.maf01.LD50_5_01.1000G_hg38.filtered.traw")))
    chr_geno_annot = chr_geno[chr_geno$SNP %in% chr_annot, ]
    rownames(chr_geno_annot) <- chr_geno_annot$SNP
    chr_geno_annot_clean <- chr_geno_annot[, 7:dim(chr_geno_annot)[2]] #select only genotype
    geno_annot_clean=rbind(geno_annot_clean, chr_geno_annot_clean)
	rm(chr_geno_annot)
  }
  print(paste0("Total num of annotated SNPs is ",dim(geno_annot_clean)[1]))
  geno_annot_r0a1 <- 2 - geno_annot_clean  ##Chek the coding for ref and alt from *.traw and *.tped file
  tmp = sub("_(.*)", "", colnames(geno_annot_r0a1))
  colnames(geno_annot_r0a1) <- tmp
  sample_size= dim(geno_annot_r0a1)[2]
  y_geno <-colSums(beta_1 * geno_annot_r0a1)
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
colnames(sim_y_df_t_IID) <- c("IID", sim_track_name)
fwrite(sim_y_df_t_IID,paste0("./simulation_values_res/sim_y_values/beta1_",beta_1,"_h2_",h2,".sim_y.txt"), sep="\t", row.names = FALSE, col.names = TRUE)
