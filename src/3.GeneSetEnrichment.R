library("readxl")

ST5=read_xlsx("F:/Dropbox/GuoLab/COREAD_ToDrGuo/CRC_TF_TWAS_Zhishan/NCrevision/Supplementary Tables_20250206.xlsx", sheet="ST5", skip=2)
ST8=read_xlsx("F:/Dropbox/GuoLab/COREAD_ToDrGuo/CRC_TF_TWAS_Zhishan/NCrevision/Supplementary Tables_20250206.xlsx", sheet="ST8", skip=2)
ST10=read_xlsx("F:/Dropbox/GuoLab/COREAD_ToDrGuo/CRC_TF_TWAS_Zhishan/NCrevision/Supplementary Tables_20250206.xlsx", sheet="ST10", skip=2)

twpsets = intersect(ST5$Genes, ST8$Genes)
threesets = intersect(twpsets, ST10$Genes)
print(threesets)

library(data.table)
library("remotes")
library("readxl")
library("enrichR")
listEnrichrSites()
setEnrichrSite("Enrichr")

dbs <- listEnrichrDbs()

#94 genes
df <- read.csv("D:/PostDoc/Projects/IBD_MBD/IBD_GWAS_TWAS/Manuscripts/DrugGenes/20250501/phase2Plus_94protein_386drugs.csv")
input <- unique(df$Genes)
print(length(input))
enriched <- enrichr(input, dbs$libraryName)

header=c("Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value", "Old.Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes")
# Initialize empty list to collect filtered data frames
filtered_list <- list()
# Loop through each enrichment result
for (sheet_name in names(enriched)) {
  df <- enriched[[sheet_name]]
  
  # Check if required column exists
  if ("Adjusted.P.value" %in% colnames(df)) {
    df_sig <- df[df$Adjusted.P.value < 0.05, ]
    
    # Proceed only if non-empty
    if (nrow(df_sig) > 0) {
      df_sig$Source <- sheet_name
      filtered_list[[sheet_name]] <- df_sig
    }
  }
}

# Combine all filtered data frames into one
if (length(filtered_list) > 0) {
  combined_df <- do.call(rbind, filtered_list)
  # Reorder columns if needed
  combined_df <- combined_df[, c(header, "Source")]
  # Save to Excel
  fwrite(combined_df, 
         file = "D:/PostDoc/Projects/IBD_MBD/IBD_GWAS_TWAS/Manuscripts/DrugGenes/20250501/enriched_genes.csv")
} else {
  message("No significant results (Adjusted.P.value < 0.05) found in any enrichment set.")
}