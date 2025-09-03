library(data.table)
library(dplyr)

#set parameters
beta1="1"
snps="1000"
wk="/CRC_TF_TWAS_zhishan2024/NC_revision_OceanCode/v2_NC_revision/"
setwd(wk)

#load sim_y files
sim_y=fread(paste0(wk, "/simulation_values_res/sim_y_values/beta1_",beta1,"_h2_0.5_SNPs_",snps,".sim_y.txt"))
sim_y_names=colnames(sim_y)[1:length(colnames(sim_y))]

#load glm res
glm_res=as.data.frame(fread(paste0(wk,"/simulation_values_res/lmer_res/beta1_",beta1,"_h2_0.5_SNPs_",snps,".sim_y.lmer.csv")))
glm_res = glm_res[!is.na(glm_res$CistromeID), ]
glm_res$CistromeID <- as.character(glm_res$CistromeID)

#load ldsc res
format_ldsc <- function(beta1) {
  ldsc_fd=paste0(wk,"/simulation_values_res/ldsc_res/beta1_", beta1,"_h2_0.5_SNPs_",snps)
  files=dir(ldsc_fd)
  ldsc_res=matrix(integer(0), nrow=length(files), ncol=6) %>% as.data.frame()
  names(ldsc_res) <- c("Beta1", "CistromeID", "Coeff", "CoeffSE", "Zvalue", "Pvalue")
  count=1
  for(file in files){
    sim_y_track= sub("\\..*$", "", file) 
    lines<-readLines(paste0(wk,"/simulation_values_res/ldsc_res/beta1_",beta1,"_h2_0.5_SNPs_",snps,"/", file))
    Cates<- Coefs <-Coefs_se <- NULL
    for(line in lines){
      if (startsWith(line, "Categories:")) {
        Cates <- strsplit(sub(".*:\\s*", "", line), "\\s+")[[1]]
      }else if (startsWith(line, "Coefficients:")) {
        Coefs <- strsplit(sub(".*:\\s*", "", line), "\\s+")[[1]]
        Coefs <- as.numeric(Coefs)
      }else if (startsWith(line, "Coefficient SE:")) {
        Coefs_se <- strsplit(sub(".*:\\s*", "", line), "\\s+")[[1]]
        Coefs_se <- as.numeric(Coefs_se)
      }
    }
    line_df = data.frame(cate=Cates,coefs=Coefs,coefs_se=Coefs_se,stringsAsFactors = FALSE)
    idx <- which(Cates == paste0(sim_y_track, "L2_0"))
    z=line_df[idx,2]/line_df[idx,3]
    p_=0
    if(z<0){ p_ <- pnorm(z)
    }else{p_ <- 1 - pnorm(z)}
    ldsc_res[count,]=c(beta1, as.character(sim_y_track), line_df[idx,2], line_df[idx,3], z, p_)
    count=count+1
  }
  return(data.frame(ldsc_res))
}

ldsc_res = format_ldsc(beta1)

out <- merge(glm_res, ldsc_res, by = "CistromeID", all.x = TRUE, suffixes = c(".glm", ".ldsc"))
out <- out[!duplicated(out), ]
fwrite(out, paste0(wk,"/simulation_values_res/lmer_ldsc_beta1_",beta1,"_h2_0.5_SNPs_",snps,".sim_y.lmer.csv"))

