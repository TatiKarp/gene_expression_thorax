library(dplyr)
library(edgeR)
library(ggplot2)
library(tibble)


setwd('/home/tatiana/Study/RP2/ATLANTIS/Umi_dedup/Dif_expr')

#data Without UMI
expression.data <- read.csv('/home/tatiana/Study/RP2/ATLANTIS/Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE)%>%
  dplyr::select(!X)

master.Table <- read.csv('./patient_data_with_asthmastatus.csv',header = TRUE)%>%
  as.data.frame() %>%
  #master.Table[!(master.Table$GenomeScan_ID %in% Exclude_list)]
  dplyr::mutate(
    GenomeScan_ID = stringr::str_trim(GenomeScan_ID), #remove whitespace from start/end of the string
    GenomeScan_ID = paste0("X", GenomeScan_ID),
    GenomeScan_ID = gsub("-", ".", GenomeScan_ID),
    age = as.numeric(age),
    sample.id = as.character(sample.id),
    smoking.status = forcats::fct_recode(
      smoking.status,
      Ex.smoker = "Ex-smoker",
      Current.smoker = "Current Smoker",
      Non.smoker= "Non-smoker"))

master.Table <- subset(master.Table, !(GenomeScan_ID %in% Exclude_list))


#Functions
select.columns.in.order <- function(dataframe, columns) {
  dataframe[, columns]
}

# Differential Expression A vs H  UMI_dedup
patients <- master.Table %>%
  dplyr::filter(
    !is.na(GenomeScan_ID)
  ) %>%
  dplyr::filter(
    !is.na(gender)
  )
sample.order <- patients %>%
  dplyr::pull(GenomeScan_ID)

expression.data <- expression.data.reduced %>%
  tibble::column_to_rownames("Gene") %>%
  select.columns.in.order(sample.order) %>%
  as.matrix()

design <- model.matrix(~0 + asthma.status + age + gender + smoking.status, data = patients)


DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL,design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]

DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")

DGEL <- edgeR::estimateDisp(DGEL, design)
#Maximizes the negative binomial likelihood to give the estimate of the common, trended and tagwise dispersions across all tags.
fit <- edgeR::glmQLFit(DGEL,design)
#Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
#Conduct genewise statistical tests for a given coefficient or contrast.

contrasts <- limma::makeContrasts(
  asthma.status = asthma.statusA - asthma.statusH, #check H-A and A-H on plot
  levels = design
)
qlf   <- edgeR::glmQLFTest(fit, contrast = contrasts[,"asthma.status"])
#Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
#Conduct genewise statistical tests for a given coefficient or contrast.

###Add gene names
library(biomaRt)
getGenedataByEnsemblId38 <- function(ensemblIds, file.location) {
  file.name <-"/home/tatiana/Study/RP2/ATLANTIS/genes_info_hg38.csv"
  if (!file.exists(file.name)) {
    if (!("mart" %in% ls())) {
      assign("mart", useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "hsapiens_gene_ensembl"
      ))
    }
    gene.list <- getBM(
      filters = "ensembl_gene_id",
      attributes = c(
        "hgnc_symbol",
        "ensembl_gene_id",
        "ensembl_transcript_id",
        "chromosome_name",
        "start_position",
        "end_position",
        "strand",
        "transcription_start_site",
        "transcript_start",
        "transcript_end",
        "external_gene_name"
      ),
      values = as.character(ensemblIds),
      mart = mart
    )
    readr::write_csv(gene.list, path = file.name)
  }
  return(
    readr::read_csv(
      file.name,
      col_types = readr::cols()
    )
  )
}

gene.data <- getGenedataByEnsemblId38(
  ensemblIds = expression.data$Gene,
  file.location = expression.dir
) %>%
  dplyr::group_by(hgnc_symbol) %>%
  dplyr::filter(
    dplyr::row_number() == 1,
    !is.na(hgnc_symbol),
    hgnc_symbol != ""
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    hgnc_symbol,
    ensembl_gene_id,
    chromosome_name,
    transcript_start,
    transcript_end
  )

#create the result table
de.results <- edgeR::topTags(qlf,  n=nrow(DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::left_join(
    y = gene.data,
    by = c("Gene" = "ensembl_gene_id")
  )# %>%
#readr::write_csv("DE.genes.ATLANTIS.csv")

#summarize:
summary(decideTests(qlf)) ##For summary use FDR<0.05!

#1*asthma.statusA -1*asthma.statusH
#Down                                   59
#NotSig                              17678
#Up                                     67

##PLOTTING
plotMD(qlf)
abline(h=c(-1,1), col="blue")


## Vulcano
library(ggrepel)

ggplot(data=de.results, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color= ifelse(FDR < 0.05,'firebrick', 'gray'))) + 
  #geom_point(aes(color= ifelse((FDR < 0.05) & (abs(logFC)>1),'firebrick', 'gray'))) + 
  theme_bw()+
  #geom_vline(xintercept=c(-1, 1), col="red") + #which log FC to choose??
  #geom_hline(yintercept=de.results$FDR <= 0.05), col="red") + #look at table
  geom_text_repel(aes(label=ifelse((FDR < 0.05 )& (abs(logFC)>1), (ifelse((!is.na(hgnc_symbol)), hgnc_symbol, Gene)), ''),
                      hjust= 0.4, vjust= 0.3))+# segment.color = 'grey50'))
  #geom_text(label=ifelse(de.results$FDR <= 0.05 & (abs(de.results$logFC)>1), de.results$hgnc_symbol , ''))+
  labs(title = "DE genes (UMI dedup)")+
  scale_color_identity(name = '', breaks= c('firebrick', 'gray') ,
                       labels = c('FDR<0.05', 'non-sign'), guide = "legend")
  xlim(-2, 5.2)


##### heatmap ######
library(pheatmap)
max.genes.in.heatmap <- 126
fc.cutoff <-0


genes.to.plot <- de.results %>%
  dplyr::select(-c(chromosome_name,transcript_start,transcript_end))%>%
  dplyr::filter( # Cutoffs
    FDR < 0.05
    #abs(logFC) > log2(fc.cutoff) ##-Inf?? Why to use this step?
  ) %>%
  dplyr::arrange(
    PValue
  ) %>%
  dplyr::filter(
    dplyr::row_number() <= max.genes.in.heatmap # Let's plot the top 35
  ) %>%
  dplyr::arrange( # Visually group the top results by Fold Change, within there sort by PValue
    logFC > 0, 
    PValue
  ) %>%
  dplyr::mutate(
    gene.change = ifelse(
      logFC > 0,
      "Up",
      "Down"))
genes.to.plot <- genes.to.plot%>%
  dplyr::mutate(
    name.to.plot = ifelse(
      FDR < 0.01,
      genes.to.plot$hgnc_symbol,
      ''))


master.Table <- master.Table %>%
  dplyr::arrange(asthma.status)

expression.data.plot.J <- expression.data.reduced %>%
  select.columns.in.order(c("Gene", master.Table$GenomeScan_ID)) #Order samples as in master.Table, since they are rearranged 
#############   Normalize with cpm  ############### 

library(matrixStats)
logcpm <- cpm(DGEL, log=TRUE)
heatmap_log_data <- as.matrix(logcpm[genes.to.plot$Gene, master.Table$GenomeScan_ID])
#color breaks
symmertic_breaks <- seq(from=-3, to= 3, length.out = 257)

col_annotation <- data.frame(Sample = master.Table$GenomeScan_ID ,Asthma.status= master.Table$asthma.status)%>%
  column_to_rownames ('Sample')
row_annotation <- data.frame(Gene=genes.to.plot$Gene, Change= genes.to.plot$gene.change)%>%
  column_to_rownames ('Gene')
jpeg(file="pheatmap_logcpm_normalization.jpg", width = 841, height = 637, res = 120)
map <- pheatmap::pheatmap(heatmap_log_data, 
                          scale='row',
                          cluster_rows = T,
                          cluster_cols = F,
                          show_rownames = T,
                          show_colnames = F,
                          color=colorRampPalette(c('blue4', 'white', 'red4'))(256),
                          breaks = symmertic_breaks,
                          annotation_col = col_annotation,
                          angle_row = 45,
                          annotation_row = row_annotation,
                          labels_row = genes.to.plot$name.to.plot,
                          fontsize_row = 4,
                          main = 'Heatmap (logcpm_normalization)')
dev.off()

###cluster within asthma ###

new_order_genes <- rownames(heatmap_log_data[map$tree_row[["order"]],])
asthma_master_table <- master.Table %>%
  filter (asthma.status=='A')
asthma_matrix <- as.matrix(logcpm[new_order_genes, asthma_master_table$GenomeScan_ID])
healthy_master_table <- master.Table %>%
  filter (asthma.status=='H')
healthy_matrix <- as.matrix(logcpm[new_order_genes, healthy_master_table$GenomeScan_ID])

#heatmap only for asthma (with preordered genes)
asthma_map <- pheatmap::pheatmap(asthma_matrix, 
                                 scale='row',
                                 cluster_rows = F,
                                 cluster_cols = T,
                                 show_rownames = T,
                                 show_colnames = F,
                                 color=colorRampPalette(c('blue4', 'white', 'red4'))(256),
                                 breaks = symmertic_breaks,
                                 annotation_col = col_annotation,
                                 angle_row = 45,
                                 annotation_row = row_annotation,
                                 #labels_row = genes.to.plot$name.to.plot,
                                 fontsize_row = 4,
                                 main = 'Heatmap (logcpm_normalization) only asthma')
#save new order of samples
new_order_asthma_samples <- colnames(asthma_matrix[,asthma_map$tree_col[["order"]]])

#heatmap only healthy (with preordered genes)
healthy_map <- pheatmap::pheatmap(healthy_matrix, 
                                  scale='row',
                                  cluster_rows = F,
                                  cluster_cols = T,
                                  show_rownames = T,
                                  show_colnames = F,
                                  color=colorRampPalette(c('blue4', 'white', 'red4'))(256),
                                  breaks = symmertic_breaks,
                                  annotation_col = col_annotation,
                                  angle_row = 45,
                                  annotation_row = row_annotation,
                                  #labels_row = genes.to.plot$name.to.plot,
                                  fontsize_row = 4,
                                  main = 'Heatmap (logcpm_normalization) only asthma')
new_order_healthy_samples <- colnames(healthy_matrix[,healthy_map$tree_col[["order"]]])

#save new order of samples
new_order_all_samples <- c(new_order_asthma_samples,new_order_healthy_samples)
new_all_matrix <- as.matrix(logcpm[new_order_genes, new_order_all_samples])

#reorder genes to plot table and master table
genes_to_plot_reordered <- genes.to.plot[match(new_order_genes, genes.to.plot$Gene),]
master_table_reordered <- master.Table[match(new_order_all_samples, master.Table$GenomeScan_ID),]

col_annotation_reodered <- data.frame(Sample = master.Table$GenomeScan_ID ,
                                      Asthma.status= master_table_reordered$asthma.status)%>% 
                                      #Smoking.status = master_table_reordered$smoking.status)
  column_to_rownames ('Sample')
row_annotation_reordered <- genes_to_plot_reordered[,c(1,8)] %>%
  remove_rownames %>% 
  column_to_rownames(var="Gene")
reordered_map <- pheatmap::pheatmap(new_all_matrix, 
                                    scale='row',
                                    cluster_rows = F,
                                    cluster_cols = F,
                                    show_rownames = T,
                                    show_colnames = F,
                                    color=colorRampPalette(c('blue4', 'white', 'red4'))(256),
                                    breaks = symmertic_breaks,
                                    annotation_col = col_annotation_reodered,
                                    angle_row = 45,
                                    annotation_row = row_annotation_reordered,
                                    labels_row = genes_to_plot_reordered$name.to.plot,
                                    fontsize_row = 4,
                                    gaps_col = 361,
                                    main = 'Heatmap (logcpm_normalization) reordered',
                                    angle_col = 315)




####################   extract the clusters   ################
###Look at the dendrogram (rows)
plot(map$tree_row)

genes_clusters_vector <- sort(cutree(map$tree_row, k=6))
genes_clusters <- data.frame(keyName=names(genes_clusters_vector), value=genes_clusters_vector, row.names=NULL)

genes_to_plot_reordered<- genes_to_plot_reordered%>%
  dplyr::left_join(
    y = genes_clusters,
    by = c("Gene" = "keyName")
  )# %>%

row_annotation_reordered <- data.frame(Gene=genes_to_plot_reordered$Gene, Change= genes_to_plot_reordered$gene.change, 
                                       Cluster=genes_to_plot_reordered$value.x)%>%
  column_to_rownames ('Gene')
# map <- pheatmap::pheatmap(heatmap_log_data, 
#                           scale='row',
#                           cluster_rows = T,
#                           cluster_cols = F,
#                           show_rownames = T,
#                           show_colnames = F,
#                           color=colorRampPalette(c('blue4', 'white', 'red4'))(256),
#                           breaks = symmertic_breaks,
#                           annotation_col = col_annotation,
#                           angle_row = 45,
#                           annotation_row = row_annotation,
#                           labels_row = genes.to.plot$name.to.plot,
#                           fontsize_row = 4,
#                           main = 'Heatmap (logcpm_normalization)')
# dev.off()

#reorder genes and samples with clusters:

reordered_map <- pheatmap::pheatmap(new_all_matrix, 
                                    scale='row',
                                    cluster_rows = F,
                                    cluster_cols = F,
                                    show_rownames = T,
                                    show_colnames = F,
                                    color=colorRampPalette(c('blue4', 'white', 'red4'))(256),
                                    breaks = symmertic_breaks,
                                    annotation_col = col_annotation_reodered,
                                    angle_row = 45,
                                    annotation_row = row_annotation_reordered,
                                    labels_row = genes_to_plot_reordered$name.to.plot,
                                    fontsize_row = 4,
                                    gaps_col = 361,
                                    gaps_row = c(18,25,59,76,86),
                                    main = 'Heatmap (logcpm_normalization) reordered')
dev.off()



write.csv(genes_clusters,'Genes_clusters_from_heatmap.csv')



