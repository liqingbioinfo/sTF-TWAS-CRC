# conda activate tftwas
library("methods")
library("glmnet")
library("stringr")
library("data.table")
library("dplyr")


"%&%" <- function(a,b) paste(a,b, sep = "")
prefix <- "model"

###### I. paramter ###### 
args <- commandArgs(trailingOnly=T)
chrom <- 22

###### II. input files ######
# snp annotation
snp_annot_file <- paste0("../data/TFs_variants_chr",chrom,".txt")
# genotype file
genotype_file <- paste0("../data/example.gt.gt") #Through dbGap
# vcf head
vcfhead <- read.table("../data/vcf.head",header = F, stringsAsFactors = F)
vcfhead <- as.character(vcfhead[1,])

# gene annotation
gene_annot_file <- "../data/gencode.v26.annotation.gene.gtf" #download from gencode URL: https://www.gencodegenes.org/human/
# gene expression
expression_file <- "../data/example_expression_matrix_rmCovPEER.tsv"


###### III. set function  ######
snp_annot <- read.table(snp_annot_file, header = T, stringsAsFactors = F)
#snp_annot$varID <- gsub("chr(\\S+)","\\1",snp_annot$varID,perl=TRUE)

get_gene_annotation <- function(gene_annot_file_name, chrom)
{
    gene_df <- read.table(gene_annot_file,header=F,stringsAsFactors =F,sep="\t",fill=T)
    gene_df1 <- filter(gene_df,V3 %in% "gene")
    geneid <- str_extract(gene_df1[,9], "ENSG\\d+.\\d+")
    genename <- gsub("gene_name (\\S+);","\\1",str_extract(gene_df1[,9], "gene_name (\\S+);"), perl=T)
    gene_used <- as.data.frame(cbind(geneid,genename,gene_df1[,c(1,4,5,3)]))
    colnames(gene_used) <- c("geneid","genename","chr","start","end","anno")
    gtf_used <- filter(gene_used,gene_used[,3] %in% ('chr' %&% chrom))
    gtf_used
}


get_gene_expression <- function(gene_expression_file_name, gene_annot) {
  expr_df <- as.data.frame(read.table(gene_expression_file_name, header = T, stringsAsFactors = F, row.names = 1))
  expr_df <- expr_df %>% select(one_of(intersect(gene_annot$geneid, colnames(expr_df))))
  expr_df <- expr_df[order(row.names(expr_df)), ]
  expr_df
}

get_maf_filtered_genotype <- function(genotype_file_name, vcfhead, maf, samples) {
  gt_df<- fread(genotype_file,sep="\t",header=F)
  colnames(gt_df) <- vcfhead
  gt_df <- as.data.frame(gt_df)
  row.names(gt_df) <- gt_df$ID
  gt_df <- gt_df[,-c(1:5)]
  gt_df1 <- as.data.frame(t(gt_df))
  gt_df1 <- gt_df1[row.names(gt_df1) %in% samples,]
  effect_allele_freqs <- colMeans(gt_df1) / 2
  gt_df1 <- gt_df1[,which((effect_allele_freqs >= maf) & (effect_allele_freqs <= 1-maf))]
#  colnames(gt_df1) <- gsub("chr(\\S+)","\\1",colnames(gt_df1),perl=TRUE)
  gt_df1 <- gt_df1[order(row.names(gt_df1)), ]
  gt_df1
}


get_gene_coords <- function(gene_annot, gene) {
  row <- gene_annot[which(gene_annot$geneid == gene),]
  c(row$start, row$end)
}


get_cis_genotype <- function(gt_df, snp_annot, coords, cis_window) {
  snp_info <- snp_annot %>% filter((pos >= (coords[1] - cis_window)  & (pos <= (coords[2] + cis_window))))
  cis_gt <- gt_df[, colnames(gt_df) %in% snp_info$varID]
  cis_gt
}

do_elastic_net <- function(cis_gt, expr_adj, n_folds, cv_fold_ids, n_times, alpha) {
    cis_gt <- as.matrix(cis_gt)
    fit <- cv.glmnet(cis_gt, expr_adj, nfolds = n_folds, alpha = alpha, keep = TRUE, type.measure='mse', foldid = cv_fold_ids[,1], parallel = FALSE)
    lambda_seq <- fit$lambda
    cvms <- matrix(nrow=length(lambda_seq), ncol=n_times)
    fits <- list()
    fits[[1]] <- fit
    cvms <- matrix(nrow = 100, ncol = n_times)
    cvms[1:length(fit$cvm),1] <- fit$cvm
    for (i in 2:(n_times)) {
      fit <- cv.glmnet(cis_gt, expr_adj, lambda = lambda_seq, nfolds = n_folds, alpha = alpha, keep = FALSE, foldid = cv_fold_ids[,i], parallel = FALSE)
      fits[[i]] <- fit
      cvms[1:length(fit$cvm),i] <- fit$cvm
    }
    avg_cvm <- rowMeans(cvms)
    best_lam_ind <- which.min(avg_cvm)
    best_lambda <- lambda_seq[best_lam_ind]
    out <- list(cv_fit = fits[[1]], min_avg_cvm = min(avg_cvm, na.rm = T), best_lam_ind = best_lam_ind, best_lambda = best_lambda)
    out
}

evaluate_performance <- function(cis_gt, expr_adj, fit, best_lam_ind, best_lambda, cv_fold_ids, n_folds) {
  n_nonzero <- fit$nzero[best_lam_ind]
  if (n_nonzero > 0) {
    R2 <- rep(0, n_folds)
    for (j in (1:n_folds)) {
      fold_idxs <- which(cv_fold_ids[,1] == j)
      tss <- sum(expr_adj[fold_idxs]**2)
      rss <- sum((expr_adj[fold_idxs] - fit$fit.preval[fold_idxs, best_lam_ind])**2)
      R2[j] <- 1 - (rss/tss)
    }
    best_fit <- fit$glmnet.fit
    expr_adj_pred <- predict(best_fit, as.matrix(cis_gt), s = best_lambda)
    tss <- sum(expr_adj**2)
    rss <- sum((expr_adj - expr_adj_pred)**2)
    
    n_samp <- length(expr_adj)
    weights <- best_fit$beta[which(best_fit$beta[,best_lam_ind] != 0), best_lam_ind]
    weighted_snps <- names(best_fit$beta[,best_lam_ind])[which(best_fit$beta[,best_lam_ind] != 0)]
    R2_mean <- mean(R2)
    R2_sd <- sd(R2)
    inR2 <- 1 - (rss/tss)
    # Old way
    pred_perf <- summary(lm(expr_adj ~ fit$fit.preval[,best_lam_ind]))
    pred_perf_rsq <- pred_perf$r.squared
    
    pred_perf_pval <- pred_perf$coef[2,4]
    out <- list(weights = weights, n_weights = n_nonzero, weighted_snps = weighted_snps, R2_mean = R2_mean, R2_sd = R2_sd,
                inR2 = inR2, pred_perf_rsq = pred_perf_rsq, pred_perf_pval = pred_perf_pval)
  } else {
    out <- list(weights = NA, n_weights = n_nonzero, weighted_snps = NA, R2_mean = NA, R2_sd = NA,
                inR2 = NA, pred_perf_rsq = NA, pred_perf_pval = NA)
  }
  out
}

do_covariance <- function(gene_id, cis_gt, rsids, varIDs, out_file) {
	model_gt <- cis_gt[,varIDs, drop=FALSE]
  	geno_cov <- cov(model_gt)
  	cov_df <- data.frame(gene=character(),rsid1=character(),rsid2=character(), covariance=double())
	for (i in 1:length(rsids)) {
	for (j in i:length(rsids)) {
	cov_df <- tryCatch(rbind(cov_df, data.frame(gene=gene_id,rsid1=rsids[i], rsid2=rsids[j], covariance=geno_cov[i,j])),
			   error = function(cond) browser())
		}
	}
	write.table(cov_df, file = out_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
}


###### IV. run analysis ######
    # parameter for analysis 
    maf=0.05
    n_times=3
    n_k_folds=10
    cis_window=1000000
    alpha=0.5

    # inite analysis
    gene_annot <- get_gene_annotation(gene_annot_file, chrom)
    expr_df <- get_gene_expression(expression_file, gene_annot)
    samples <- rownames(expr_df)
    n_samples <- length(samples)
    genes <- colnames(expr_df)
    n_genes <- length(expr_df)
    gt_df <- get_maf_filtered_genotype(genotype_file, vcfhead, maf, samples)
 
    samples <- intersect(row.names(expr_df),row.names(gt_df))
    expr_df <- expr_df[row.names(expr_df) %in% samples,]
    #dim(expr_df)
    n_samples <- length(samples)

    seed <- NA
    seed <- ifelse(is.na(seed), sample(1:1000000, 1), seed)
    set.seed(seed)
    cv_fold_ids <- matrix(nrow = n_samples, ncol = n_times)
    for (j in 1:n_times)
    cv_fold_ids[,j] <- sample(1:n_k_folds, n_samples, replace = TRUE)
  
    # output file
    model_summary_file <- '../twasres/' %&% prefix %&% '_chr' %&% chrom %&% '_model_summaries.txt'
    model_summary_cols <- c('gene_id', 'gene_name', 'alpha', 'cv_mse', 'lambda_iteration', 'lambda_min',
                            'n_snps_in_model','cv_R2_avg', 'cv_R2_sd', 'in_sample_R2', 'pred_perf_R2',
                            'pred_perf_pval')
    write(model_summary_cols, file = model_summary_file, ncol = 12, sep = '\t')
  
    weights_file <- '../twasres/' %&% prefix %&% '_chr' %&% chrom %&% '_weights.txt'
    weights_col <- c('gene_id', 'rsid', 'varID', 'ref', 'alt', 'beta')
    write(weights_col, file = weights_file, ncol = 6, sep = '\t')
  
    covariance_file <- '../twasres/' %&% prefix %&% '_chr' %&% chrom %&% '_covariances.txt'
    covariance_col <- c('GENE', 'RSID1', 'RSID2', 'VALUE')
    write(covariance_col, file = covariance_file, ncol = 4, sep = ' ')
  
  for (i in 1:n_genes) {

    cat(i, "/", n_genes, "\n")
    gene <- genes[i]
    gene_name <- as.character(gene_annot$genename[gene_annot$geneid == gene])
    model_summary <- c(gene, gene_name, alpha, NA, NA, NA, 0, NA, NA, NA, NA, NA)
    coords <- get_gene_coords(gene_annot, gene)
    cis_gt <- get_cis_genotype(gt_df, snp_annot, coords, cis_window)

    tryCatch({

    if (ncol(cis_gt) >= 2) {
      adj_expression <- expr_df[,i]
      elnet_out <- do_elastic_net(cis_gt, adj_expression, n_k_folds, cv_fold_ids, n_times, alpha)
      if (length(elnet_out) > 0) {
        eval <- evaluate_performance(cis_gt, adj_expression, elnet_out$cv_fit, elnet_out$best_lam_ind, elnet_out$best_lambda, cv_fold_ids, n_k_folds)
        model_summary <- c(gene, as.character(gene_name), alpha, elnet_out$min_avg_cvm, elnet_out$best_lam_ind,
                           elnet_out$best_lambda, eval$n_weights, eval$R2_mean, eval$R2_sd, eval$inR2,
                           eval$pred_perf_rsq, eval$pred_perf_pval)
        if (eval$n_weights > 0) {
          weighted_snps_info <- snp_annot %>% filter(varID %in% eval$weighted_snps) %>% select(SNP,varID, ref , effect)
          if (nrow(weighted_snps_info) == 0)
            browser()
          weighted_snps_info$gene <- gene
          weighted_snps_info <- weighted_snps_info %>% merge(data.frame(weights = eval$weights, varID=eval$weighted_snps), by = 'varID') %>% select(gene, SNP, varID, ref, effect, weights)
          write.table(weighted_snps_info, file = weights_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
	  do_covariance(gene, cis_gt, weighted_snps_info$SNP, weighted_snps_info$varID, covariance_file)
        }
      }
    }
  write(model_summary, file = model_summary_file, append = TRUE, ncol = 12, sep = '\t')

    },error=function(e){})
    
  }
