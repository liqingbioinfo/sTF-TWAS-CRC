#conda activate tftwas
library("sqldf")
library("stringr")
library("data.table")
library("dplyr")
args <- commandArgs(trailingOnly=T)


###### I. paramter ###### 
model_prefix <- "model"
wk="../results/"


######II. Collect all model summary files [ensure one ENSG name]
model_summary <- matrix(integer(0), nrow = 0, ncol = 6) %>% as.data.frame()
names(model_summary) <- c("gene_id","genename","pred.perf.R2","n.snps.in.model","pred.perf.pval","pred_perf_pval.1")

model_weights <- matrix(integer(0), nrow = 0, ncol = 5) %>% as.data.frame()
names(model_weights) <- c("rsid","gene_id","beta","ref","alt") #gene_id rsid    varID   ref     alt     beta

model_covars <- matrix(integer(0), nrow = 0, ncol = 4) %>% as.data.frame()
names(model_covars) <- c('gene_id', 'rsid1', 'rsid2', 'corvarianceValues')

####Collect model summary#####
for(i in 1:22){
	chrom=as.character(i)
	if (file.exists(paste0(wk, model_prefix, "_chr",chrom , "_model_summaries.txt"))) {
		model_summary_perchr=as.data.frame(fread(paste0(wk, model_prefix, "_chr",chrom , "_model_summaries.txt")))
		if(dim(model_summary_perchr)[2]==12){
			model_summary_perchr_remain=model_summary_perchr[, c(1,2,11,7,12,12)] #$F[0]\t$F[1]\t$F[10]\t$F[6]\t$F[11]\t$F[11]
			model_summary = rbind(model_summary, model_summary_perchr_remain)
		}else{
			cat("Error: reading ",chrom," model summary file does not have 12 columns! Exit")
			q()
		}
	} else {
	  cat("Warning: ", model_prefix, "_chr", as.character(chrom), "_model_summaries.txt does not exit")
	}
}
##Remove pred.perf.R2 is NA and keep only one ENSG
model_summary <- model_summary[!is.na(model_summary$pred_perf_R2),]
model_summary_unique <- model_summary %>% arrange(pred_perf_R2) %>% distinct(gene_id, .keep_all=TRUE)
model_summary_unique <- model_summary_unique %>% distinct(gene_name, .keep_all=TRUE)
colnames(model_summary_unique) <- c("gene", "genename" ,"pred.perf.R2", "n.snps.in.model", "pred.perf.pval", "pred.perf.qval")
fwrite(model_summary_unique, paste0(wk, model_prefix, "_model_summaries.csv"), row.names = FALSE, col.names = TRUE)

####Collect model weights#####
for(i in 1:22){
	chrom=as.character(i)
	if (file.exists(paste0(wk, model_prefix, "_chr",chrom , "_weights.txt"))) {
		model_weights_perchr=as.data.frame(fread(paste0(wk, model_prefix, "_chr",chrom , "_weights.txt")))
		if(dim(model_weights_perchr)[2]==6){
			model_weights_perchr_remain=model_weights_perchr[, c(2,1,6,4,5)] #"$F[1],$F[0],$F[5],$F[3],$F[4]"
			model_weights = rbind(model_weights, model_weights_perchr_remain)
		}else{
			cat("Error: reading ",chrom," model weights file does not have 6 columns! Exit")
			q()
		}
	}else{
		cat("Warning: ", model_prefix, "_chr", as.character(chrom), "_weights.txt does not exit")
	}
}
model_weights <- model_weights[!is.na(model_weights$beta),]
model_weights_unique <- model_weights %>% distinct(gene_id, rsid , .keep_all=TRUE)
colnames(model_weights_unique) <- c("rsid","gene","weight","ref_allele","eff_allele")
fwrite(model_weights_unique, paste0(wk, model_prefix, "_model_weights.csv"), row.names = FALSE, col.names = TRUE)

####Collect model covariates#####
for(i in 1:22){
	chrom=as.character(i)
	if (file.exists(paste0(wk, model_prefix, "_chr",chrom , "_covariances.txt"))) {
		model_covars_perchr=as.data.frame(fread(paste0(wk, model_prefix, "_chr",chrom , "_covariances.txt")))
		if(dim(model_covars_perchr)[2]==4){
			model_covars = rbind(model_covars, model_covars_perchr)
		}else{
			cat("Error: reading ",chrom," model covariances file does not have 6 columns! Exit")
			q()
		}
	}else{
		cat("Warning: ", model_prefix, "_chr", as.character(chrom), "_covariances.txt does not exit")
	}
}
model_covars <- model_covars[!is.na(model_covars$corvarianceValues),]
model_covars_unique <- model_covars %>% distinct(gene_id, rsid1, rsid2, .keep_all=TRUE)
colnames(model_covars_unique) <-c("GENE","RSID1","RSID2","VALUE")
fwrite(model_covars_unique, paste0(wk, model_prefix, "_cov.txt"), sep=" ", row.names = FALSE, col.names = TRUE)
system(paste("gzip", paste0(wk, model_prefix, "_cov.txt")))


######III. Generate a db file for SPrediXcan
Summary <- model_summary_unique
Weight <- model_weights_unique
dbfile <- paste0(wk, model_prefix, ".db")

db <- dbConnect(SQLite(), dbname= dbfile)
dbWriteTable(conn = db, name = "extra", value = Summary, row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "weights", value = Weight,row.names = FALSE, header = TRUE)

