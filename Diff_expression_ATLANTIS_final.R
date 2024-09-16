## This script will perform DGE analysis for ATLANTIS dataset comparing asthma vs healthy subjects 
library(dplyr)
library(edgeR)
library(ggplot2)
library(tibble)
library(biomaRt)

setwd("/Users/tatiana/Work/RP2/ATLANTIS")

# metafile 
master.Table <- read.csv("./Umi_dedup/Dif_expr/ATLANTIS_master_table_QC.csv", header =TRUE) %>%
  filter (QC_check == 'YES')

# UMI deduplicated read count table
expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE) %>%
  tibble::column_to_rownames("Gene") %>%
  dplyr::select(c(master.Table$GenomeScan_ID)) %>%
  as.matrix()

# create design martix, covariates include sex, age, smoking status
design <- model.matrix(~0 + asthma.status + age + gender + smoking.status, data = master.Table)

# perform DE analysis
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL,design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
DGEL <- edgeR::estimateDisp(DGEL, design)
fit <- edgeR::glmQLFit(DGEL, design, legacy = TRUE) 
# legacy = TRUE: for edgeR 4.0 and above- to reproduce the results of edgeR 3.0)

# define contrast for the comparison
contrasts <- limma::makeContrasts(
  asthma.status = asthma.statusA - asthma.statusH,
  levels = design
)
qlf   <- edgeR::glmQLFTest(fit, contrast = contrasts[,"asthma.status"])

# Add gene names
## add gene names 
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_ids <- rownames(DGEL$counts)  # here are the ensembl gene IDs
all_new_gene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id", values = gene_ids, mart = ensembl)

#create the result table
de.results <- edgeR::topTags(
  qlf,
  n=nrow(DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::left_join(
    y = all_new_gene,
    by = c("Gene" = "ensembl_gene_id")) %>%
  mutate (logFC = round (logFC, 2),
          PValue = signif(PValue, digits = 3),
          FDR = signif(FDR, digits = 3))# %>%
#readr::write_csv("./Umi_dedup/Dif_expr/DE.genes.ATLANTIS.csv")

#summarize:
summary(decideTests(qlf)) 

#1*asthma.statusA -1*asthma.statusH
#Down                                   59
#NotSig                              17678
#Up                                     67




