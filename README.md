# sTF-TWAS-CRC

## Overview
In this work, we analyzed 218 TF chromatin immunoprecipitation sequencing (ChIP-Seq) datasets from colorectal cancer (CRC) related cells, along with GWAS data comprising 100,204 CRC cases and 154,587 controls of East Asian and European ancestries, using the generalized mixed models we developed. Under our sTF-TWAS framework, we performed trans-ancestry analyses to evaluate associations between genetically predicted gene expression, alternative splicing (AS), and alternative polyadenylation (APA) with CRC risk, using RNA-seq data generated in normal colon tissue samples obtained from 364 individuals of Asian-ancestry and 707 individuals of European-ancestry. Differential gene expression analyses across colon normal, adenoma and carcinoma tissues were conducted using single cell RNA-seq data from over 100 individuals. Gene functional validation was performed using CRISPR-Cas9 knockdown in CRC cell lines.

![My Image](./Figures/Supplementary_Figure_1.png)

**Step1:** Identify the risk TFs for CRC.

**Step2:** Constructions of sTFTWAS models using gene expression, alternative splicing and Alternative polyadenylation (APA) for European ancestry and Asian ancestry.

**Step3** Associations sTFTWAS models with colorectal cancer GWAS summary statistics from large population to identify risk proteins for CRC.

**Step4** Differential expression analysis from single cell data

## Methods
### Identification of CRC risk TFs
We used our recently developed generalized mixed model approach7 to investigate associations of CRC risk with variations of TF-DNA binding by a single TF. Please refer to our previous publications [PMID: 34518541](https://pubmed.ncbi.nlm.nih.gov/34518541/) and github repository(https://github.com/XingyiGuo/BC-TFvariants/) for more information.

### Constructions of sTFTWAS models
We utilized the weight matrix and the summary statistics from CRC GWAS datasets consisting of 186,072 individuals of European ancestry and 72,272 individuals of East Asian ancestry, we evaluated the association between gene expression (or AS, APA) and CRC risk under the sTF-TWAS framework.

- Executive code: xxx
- Input files:
1)
2)
- Output files

### Identification of CRC risk events through associations
We utilized the weight matrix and the summary statistics from CRC GWAS datasets consisting of 186,072 individuals of European ancestry and 72,272 individuals of East Asian ancestry, we evaluated the association between gene expression (or AS, APA) and CRC risk under the sTF-TWAS framework.

- Executive code: xxx
- Input files:
1)
2)
- Output files

### Identification of CRC risk events through associations
We conducted differential expression analysis across colorectal cancer (CRC) developmental stages to identify key genes associated with tumorigenesis. Using a pseudobulk gene count matrix derived from single-cell RNA sequencing (scRNA-seq) data, we applied DESeq2 to detect differentially expressed genes between normal tissues and precancerous/cancerous stages along both the normal-serrated polyp-carcinoma and normal-adenoma-carcinoma pathways.

- Executive code: xxx
- Input files:
1)
2)
- Output files

## Contact
Qing Li: qing.li@vumc.org
Xingyi Guo: xingyi.guo@vumc.org
Zhishan Chen: zhishan118@gmail.com
