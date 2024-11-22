##### TCGA analysis in this paper - R Source code #####

#Set seed for reproducibility
set.seed(8)

#Load necessary packages
library(TCGAbiolinks) 
library(tidyverse)
library(SummarizedExperiment) 
library(CIBERSORT)
library(IOBR)
library(ggExtra)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readxl)
library(data.table) 
library(ggsignif)
library(limma)

###Figure S1B###
setwd("########")
project <- "TCGA-SKCM" 
data_category <- "Transcriptome Profiling"
data_type <- "Gene Expression Quantification"
workflow_type <- "STAR - Counts"
query <- GDCquery(project = project,
                  data.category = data_category,
                  data.type = data_type,
                  workflow.type = workflow_type)
GDCdownload(query = query,files.per.chunk = 50)
GDCprepare(query,save = T,save.filename = paste0(project,"_transcriptome.Rdata"))
# load data
load(file = paste0(project,"_transcriptome.Rdata"))
expr_counts_mrna <- assay(se_mrna,"unstranded")
expr_tpm_mrna <- assay(se_mrna,"tpm_unstrand")
expr_fpkm_mrna <- assay(se_mrna,"fpkm_unstrand")
expr_fpkm_lnc <- assay(se_lnc,"fpkm_unstrand")
mrna <- assay(se_mrna,"fpkm_unstrand")   
symbol_mrna <- rowData(data)$gene_name
expr_counts_mrna_symbol <- cbind(data.frame(symbol_mrna),
                                  as.data.frame(mrna))
expr_read <- expr_counts_mrna_symbol %>% 
  as_tibble() %>% 
  mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
  arrange(desc(meanrow)) %>% 
  distinct(symbol_mrna,.keep_all=T) %>% 
  select(-meanrow) %>% 
  column_to_rownames(var = "symbol_mrna") %>% 
  as.data.frame()

cibersort <- deconvo_tme(eset = expr_read, 
                           method = "cibersort",
                           arrays = F,
                           perm = 1000)
cibersort <- counts %>% column_to_rownames("ID")
genes <- c("HPGDS")
genes_expr <- as.data.frame(t(expr_read[rownames(expr_read) %in% genes,]))
colnames(genes_expr)[1] <- 'ID'
df2 <- merge(genes_expr, cibersort1, by = "row.names", all = TRUE)
# plot
p12 <- ggplot(df2, aes(HPGDS,df2$T_cells_CD8_CIBERSORT))+
     geom_point(col="black")+
     geom_smooth(method=lm,se=T, na.rm=T, fullrange=T, size=1, col="#006fbc")+
     theme_minimal()+
     xlab("HPGDS")+
     ylab("CD8 T cell")+
     labs(title="TCGA-SKCM")+
     theme(legend.position = "", 
           legend.title = element_text(face="bold", size=25),
           legend.text = element_text(size=10),
           axis.text.y = element_text(size = 10),
           axis.text.x = element_text(size = 10, angle = 0),        
           axis.title = element_text(face="bold", size=20),
           axis.title.x = element_text(color = "black", size = 20, face = "italic"),
           plot.title = element_text(size = 25, face ="bold", hjust = 0.5,
                                     margin = margin(t = 10, b = 20, unit ="pt")))+
     stat_cor(method = "spearman", color = "black",label.x = 6, label.y = 0.7,size=5)
#add marginal histogram
p13 <- ggMarginal(p12, type = "histogram", xparams = list(fill ="orange"),
            yparams = list(fill ="skyblue"))

ggsave(p13, file="macrophage.jpg", height=8, width=8)


###Figure S11F###
setwd("########")
PAAD.fpkm <- fread("TCGA-PAAD.htseq_fpkm.tsv.gz", header = T, sep = '\t', data.table = F)
PAAD.pro <- fread("gencode.v22.annotation.gene.probeMap", header = T, sep = '\t', data.table = F)
PAAD.pro <- PAAD.pro[ , c(1, 2)]
PAAD.fpkm.pro <- merge(PAAD.pro, PAAD.fpkm, by.y  = "Ensembl_ID", by.x = "id" )
PAAD.fpkm.pro <- distinct(PAAD.fpkm.pro,gene, .keep_all = T)
rownames(PAAD.fpkm.pro) <- PAAD.fpkm.pro$gene
PAAD.fpkm.pro <- PAAD.fpkm.pro[ , -c(1,2)]

PAAD.phe <- fread("TCGA-PAAD.GDC_phenotype.tsv.gz", header = T, sep = '\t', data.table = F)
PAAD.phe$submitter_id.samples[1:5]
rownames(PAAD.phe) <- PAAD.phe$submitter_id.samples
table(PAAD.phe$sample_type.samples)
PAAD.phe.t <- filter(PAAD.phe, sample_type.samples == "Primary Tumor")
#extract primary tumor and normal samples
PAAD.phe.n <- filter(PAAD.phe, sample_type.samples == "Solid Tissue Normal")
#merge primary tumor and normal samples
PAAD.phe.t <- rbind(PAAD.phe.t, PAAD.phe.n)
merge_phe_fpkm <- intersect(rownames(PAAD.phe.t), colnames(PAAD.fpkm.pro))
PAAD.exp <- PAAD.fpkm.pro[ , merge_phe_fpkm]
#GETX data
gtex.exp <- fread("gtex_RSEM_gene_fpkm.gz", header = T, sep = '\t', data.table = F)
gtex.pro <- fread("gencode.v22.annotation.gene.probemap", header = T, sep = '\t', data.table = F)
gtex.pro <- gtex.pro[, c(2,3)]
gtex.fpkm.pro <- merge(gtex.pro, gtex.exp, by.y ="sample", by.x = "id" )
gtex.phe <- fread("GTEX_phenotype.gz", header = T, sep = '\t', data.table = F)
rownames(gtex.phe) <- gtex.phe$Sample
colnames(gtex.phe) <- c("Sample", "body_site_detail (SMTSD)", "primary_site", "patient", "cohort")
#take pancreas
gtex.phe.s <- filter(gtex.phe, primary_site == "Pancreas")
merge_phe_fpkm_gtex <- intersect(rownames(gtex.phe.s), colnames(gtex.fpkm.pro)) 
gtex.s <- gtex.fpkm.pro[ , c("gene", merge_phe_fpkm_gtex)]
gtex.s <- distinct(gtex.s, gene, .keep_all = T)
rownames(gtex.s) <- gtex.s$gene
#gtex data is log2(fpkm+0.001)
gtex.s <- gtex.s[ , -1]
gtex.s2 <- 2^gtex.s
gtex.s3 <- log2(gtex.s2-0.001+1)
#merge data and remove batch effect
all.data <- merge(gtex.s3, PAAD.exp, by = 0)
all.data <- column_to_rownames(all.data, "Row.names")
nromalized.data <- normalizeBetweenArrays(all.data)
nromalized.data <- as.data.frame(nromalized.data)
exprSet.all.r=nromalized.data[c("HPGDS"),]
exprSet.all.r=t(exprSet.all.r)
exprSet.all.r=as.data.frame(exprSet.all.r)
write.csv(exprSet.all.r, file = "HPGDS.csv")
#Remove 4 adjcant samples in TCGA
x=c(rep("Normal",167),rep("Tumor",177))
exprSet.all.r$Type=x

p20 <- ggboxplot(exprSet.all.r, x = "Type", y = "HPGDS", 
          color = "Type", palette = c("#00AFBB", "#F8766D"),
          add = "jitter", ylab = "HPGDS expression", xlab = "")+
    stat_compare_means(comparisons = list(c("Normal", "Tumor")), label = "p.adj")+
    theme_classic()+theme(legend.position = "none")
ggsave("HPGDS.jpg", p20, width = 4, height = 4, dpi = 300)
