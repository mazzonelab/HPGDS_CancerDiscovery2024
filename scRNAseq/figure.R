##### scRNAseq in this paper - R Source code #####

#Set seed for reproducibility
set.seed(6)

#Load necessary packages
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(tidyverse)
library(scCustomize)
library(Hmisc)
library(SingleR)
library(EnhancedVolcano)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(stats)
library(purrr)

###Figure 1A###
setwd("########")
mac <- readRDS("macrophage_cluster.rds"))
Idents(mac)=mac$`BT/OT`
OT <- subset(mac, idents = "OT")
BT <- subset(mac, idents = "BT")
object.markers <- FindMarkers(OT, ident.1 = "NR", ident.2 = "R", group.by = "Response", logfc.threshold = 0, min.pct = 0, pseudocount.use=0.01)
marker <- c("CBR3","CBR1","HPGD","AKR1C3","GPX4","FAM213B","ALOX12","ALOX15","ALOX15B","ALOX5",
"CYP4F8","CYP1A1","CYP1A2","CYP1B1","CYP3A7","CYP2A13","CYP2C8","CYP2C18",
"CYP2D6","CYP2F1","CYP2J2","CYP3A5","CYP4A11","CYP4B1","CYP4Z1","CYP4X1","CYP4A22","CYP2S1","CYP4F11","CYP3A43",
"CYP4F12","CYP2U1","CYP4F2","CYP4F3","PTGS1","PTGS2","PTGR1","PTGR2","PLA2G4B","PLA2G4E","PLB1",
"PLA2G4F","PLA2G2D","PLA2G2C","PLA2G1B","PLA2G2A","PLA2G5","PLA2G2F","PLA2G12A","PLA2G6","PLA2G10",
"PLA2G12B","PLA2G4C","JMJD7-PLA2G4B","HRASLS2","PLA2G16","EPHX2","EPHX3","LTA4H","GGT1","GGT5","GGTLC2","GGTLC1","ALOXE3",
"LTC4S","HPGDS","PTGDS","PTGES3","PTGES2","PTGES","PTGIS","TBXAS1")

object <- object.markers[rownames(object.markers) %in% marker,]    ##marker in gene list)
 
cluster1.markers <- object %>%
     mutate(Difference = pct.1 - pct.2) %>% 
     rownames_to_column("gene")
 
cluster1.markers$threshold="ns"
cluster1.markers[which(cluster1.markers$avg_log2FC  >= 1.5 & cluster1.markers$p_val_adj <0.05),]$threshold="up";
cluster1.markers[which(cluster1.markers$avg_log2FC  <= -1.5 & cluster1.markers$p_val_adj < 0.05),]$threshold="down";
cluster1.markers$threshold=factor(cluster1.markers$threshold, levels=c('down','ns','up'))

p1 <- ggplot(cluster1.markers, aes(x=Difference, y=avg_log2FC, color=threshold)) + 
     geom_point(size=2) + 
     scale_color_manual(name="differential expression",
                        values=c("blue","grey","red")) + 
     geom_label_repel(data=subset(cluster1.markers, avg_log2FC >= 1.5 &  p_val_adj <= 0.05), 
                      aes(label=gene),  
                      color="black", 
                      segment.colour = "black",
                      segment.size = 0.5,  
                      size=7)+
     geom_label_repel(data=subset(cluster1.markers, avg_log2FC <= -1.5 &  p_val_adj <= 0.05), 
                      aes(label=gene), label.padding = 0.1, 
                      color="black",
                      segment.colour = "black",
                      segment.size = 0.5, size=7)+
     geom_vline(xintercept = 0.0,linetype=4)+
     geom_hline(yintercept = 0,linetype=4)+
     theme_classic()+theme(axis.line = element_line(color = "black", size = 1),
                           axis.text.y = element_text(size = 12),
                           axis.text.x = element_text(size = 12,hjust = 1,vjust = 1))
 
ggsave(p1,filename = "volcano_OT.jpg",width = 8,height = 6)


###Figure 1B###
setwd("########")
GC_immune <- readRDS("GC_immune.rds"))
p2 <- FeaturePlot_scCustom(seurat_object = GC_immune, features = "HPGDS", order = F, pt.size = 1.2, alpha_exp = 1)+
     theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
     scale_x_continuous("UMAP 1")+scale_y_continuous("UMAP 2")+
     theme_bw()+ 
     theme(legend.position = c(0.95, 0.5),
           panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
           plot.title = element_text(hjust = 0.5,size=25, face="bold"),
           axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),
           axis.title.y = element_text(color = "black", size = 20),
           axis.title.x = element_text(color = "black", size = 20)
     )
ggsave(p2, file="FeaturePlot.jpg", width = 12, height = 15)


###Figure 1C###
setwd("########")
counts <- read.csv(file = "GSE115978_tpm.csv", header = TRUE, row.names = 1)
metadata <- read.csv("GSE115978_metadata.csv")
metadata <- metadata %>% column_to_rownames("cells")
GSE115978 <- CreateSeuratObject(counts = counts, meta.data = metadata, project = "mel", min.cells = 3, min.features = 100)
Idents(GSE115978)=GSE115978$cell.types
GSE115978 <- subset(GSE115978, idents = c("Macrophage","Endo.","T.CD4","CAF","T.CD8","T.cell","NK","B.cell"))
p3 <- VlnPlot(GSE115978, features = "HPGDS", add.noise = F,cols = c("#F8766D","#CD9600","#7CAE00","#0CB702","#00B8E7","#00A9FF","#C77CFF","#FF61CC"))+scale_x_discrete(limits = c("B.cell","CAF","Endo.","Macrophage","NK","T.CD4","T.CD8", "T.cell")) +
    NoLegend()+xlab("")  
ggsave(p3, file="HPGDS_GSE1195978.jpg", width = 8, height = 6)


###Figure 1D###
setwd("########")
count=read.delim("GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz", stringsAsFactors = F, header = T)
count[1:10,1:10]
tirosh_genes <- count[-1,-1]
tirosh_genes[1:5,1:5]
meta <- read.delim("meta.txt")
meta <- meta %>% column_to_rownames("title")
GSE120575 <- CreateSeuratObject(counts = tirosh_genes, meta.data = meta, min.cells = 3, min.features = 100)
#remove NA columns
toRemove <- "H9_P5_M67_L001_T_enriched"
GSE120575 <- GSE120575[,!colnames(GSE120575) %in% toRemove]
Idents(GSE120575)=GSE120575$characteristics..response
GSE120575 <- subset(GSE120575, idents = c("Non-responder", "Responder"))
#downstream analysis
GSE120575 <- NormalizeData(GSE120575, normalization.method = "LogNormalize", scale.factor = 10000)
GSE120575 <- FindVariableFeatures(GSE120575, selection.method = "vst", nfeatures = 2000)
GSE120575 <- ScaleData(GSE120575)
GSE120575 <- RunPCA(GSE120575)
GSE120575 <- FindNeighbors(GSE120575, dims = 1:10)
GSE120575 <- FindClusters(GSE120575, resolution = 0.4)
GSE120575 <- RunUMAP(GSE120575, dims = 1:10)
#rename clusters
new.cluster.ids <- c("T_cells", "T_cells", "T_cells", 
                     "B_cell", "NK_cell", "T_cells",
                     "Monocyte", "Macrophage", "B_cell","B_cell","B_cell")
names(new.cluster.ids) <- levels(GSE120575)
GSE120575 <- RenameIdents(GSE120575, new.cluster.ids)
GSE120575$celltype <- Idents(GSE120575)
##plot vlnplot
p4 <- VlnPlot(GSE120575, features = "HPGDS", add.noise = T, cols = c("#F8766D","#CD9600","#00A9FF","#FF61CC","#F8766D"))+scale_x_discrete(limits = c("B_cell","NK_cell","Monocyte","Macrophage","T_cells")) +
     NoLegend()+xlab("")  
ggsave(p4, file="HPGDS_GSE120575.jpg", height = 8, width = 6)


###Figure 1E###
marker <- c("HPGDS","LYVE1","MAF","MRC1","SELENOP","FOLR2","NRP1","CD209","SLC40A1","CD28","FGFR1","CD163L1","CD163","TNFRSF1A","S100A4","S100A6","MARCO")
df <- DotPlot(total, features = marker)$data 
p5 <- ggplot(df, aes(x=features.plot,y = id,size = pct.exp, color = avg.exp.scaled))+
     geom_point() + 
     scale_size("Percent\nExpressed", range = c(0,8)) +
     scale_y_discrete(position = "left") +
     scale_color_gradientn(colours = brewer.pal(5, "Reds"),
                           guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                           name = "Average\nexpression") +
     cowplot::theme_cowplot() +
     ylab("") + xlab("") + 
     theme_bw(base_rect_size = 2, base_line_size = 0.5) +
     theme(
         axis.text.x = element_text(size=15, angle=45, hjust=1, color="black", lineheight = 10, face = 2),
         axis.text.y = element_text(size=12, color="black", face = 2),
         axis.title = element_text(size=14, face = 2),
     )+coord_flip()   
p6<-  p5 +scale_y_discrete(limits = c("Macrophages_LYVE1","Macrophages_PLA2G2D","Macrophages_SPP1","Macrophages_CCL3","Macrophages_CXCL9","Monocytes"))
ggsave(p6,filename = "dotplot_cluster.jpg",width = 6,height = 12)


###Figure 1F###
high <- WhichCells(mac, expression = HPGDS > 0.5)
low <- WhichCells(mac, expression = HPGDS < 0.3)
subset_cells <- c(high, low)
mac_sub <- subset(x = mac, cells = subset_cells)
mac_sub@meta.data[high, "HPGDS"] <- "higher"
mac_sub@meta.data[low, "HPGDS"] <- "lower"
diff <- FindMarkers(mac_sub, min.pct = 0.1, group.by = "HPGDS", 
                     logfc.threshold = 0.1,
                     ident.1 ="higher",
                     ident.2="lower")  
diff <- diff[-1,]
p7 <- EnhancedVolcano(diff,
                      lab = rownames(diff),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      pCutoff = 0.05,
                      FCcutoff = 0.1,
                      pointSize = 1.2,
                      labSize = 2.5,
                      boxedLabels = T,
                      drawConnectors = T,
                      widthConnectors = 0.5,
                      labFace = 'bold',
                      selectLab = c("LYVE1","MRC1","LYZ","HLA-B","HLA-A","STAT1","SOD2","CXCL9","ALOX5AP",
                                    "NRP1","CD28","MAF","FCGBP","FCN1","HLA-C","CXCL10","CD276","MSR1","MS4A4A","VSIG4","C1QB"))

ggsave(p7,filename = "volcano_DEG.jpg",width = 5,height = 7)


###Figure 1I###
marker <- c("HPGDS","LYVE1","MAF","MRC1","SELENOP","FOLR2","NRP1","CD209","SLC40A1","CD28","FGFR1","CD163L1","CD163","TNFRSF1A","S100A4","S100A6","MARCO")
mac@meta.data <- tidyr::unite(mac@meta.data,"RT", "Response", "BT/OT",sep="_",remove=FALSE)
sample_table <- as.data.frame(table(mac@meta.data$RT))
sample_table <- sample_table %>%
    arrange(desc(Var1)) %>% 
    mutate(lab.ypos = cumsum(Freq) - 0.5*Freq)

df <- DotPlot(mac, features = marker, group.by = "RT")$data  
p8<- ggplot(df, aes(x=features.plot,y = id,size = pct.exp, color = avg.exp.scaled))+
    geom_point() + 
    scale_size("Percent\nExpressed", range = c(0,8)) +
    scale_y_discrete(position = "left") +
    scale_color_gradientn(colours = brewer.pal(5, "Reds"),
                          guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                          name = "Cell number
NR_BT=575
R_BT=499
NR_OT=812
R_OT=250
    
Average\nexpression") +
    cowplot::theme_cowplot() +
    ylab("") + xlab("") + 
    theme_bw(base_rect_size = 2, base_line_size = 0.5) +
    theme(
        axis.text.x = element_text(size=15, angle=45, hjust=1, color="black", lineheight = 10, face = 2),
        axis.text.y = element_text(size=12, color="black", face = 2),
        axis.title = element_text(size=14, face = 2),
    )+coord_flip()  
 
p9<- p8+ scale_y_discrete(limits = c("NR_BT","R_BT","NR_OT","R_OT"))
ggsave(p9,filename = "dotplot_treatment111.jpg",width = 6,height = 12)


###Figure 2A###
GSE146613 <- Read10X(data.dir = "path/to/object/")
GSE146613[["percent.mt"]] <- PercentageFeatureSet(GSE146613, pattern = "^mt-")
GSE146613 <- CreateSeuratObject(counts = tumor, project = "GSE1466133k", min.cells = 3, min.features = 100)
GSE146613 <- CreateSeuratObject(counts = GSE146613, project = "GSE1466133k", min.cells = 3, min.features = 100)
GSE146613[["percent.mt"]] <- PercentageFeatureSet(GSE146613, pattern = "^mt-")
VlnPlot(GSE146613, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
GSE146613 <- subset(GSE146613, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
GSE146613 <- NormalizeData(GSE146613, normalization.method = "LogNormalize", scale.factor = 10000)
GSE146613 <- FindVariableFeatures(GSE146613, selection.method = "vst", nfeatures = 2000)
GSE146613 <- ScaleData(GSE146613)
GSE146613 <- RunPCA(GSE146613)
GSE146613 <- FindNeighbors(GSE146613, dims = 1:20)
GSE146613 <- FindClusters(GSE146613, resolution = 0.4)
GSE146613 <- RunUMAP(GSE146613, dims = 1:20)
GSE146613.markers <- FindAllMarkers(GSE146613, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_genes <- GSE146613.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
new.cluster.ids <- c("Immunostimulatory_Mac", "Neutrophils", "Monocytes", "Immunosuppressive_Mac", "pDC2",
"CD4_Tcells", "NK", "Macrophages_Top2a", "CD8_Tcells",
"Migratory_DCs", "moDCs", "pDCs", "Stromal","Mast_cells","Skin_resident_Macrophages","B_cells","Neutrophils_Camp")
GSE146613 <- RenameIdents(GSE146613, new.cluster.ids)
names(new.cluster.ids) <- levels(GSE146613)
GSE146613$celltype <- Idents(GSE146613)
view(GSE146613@meta.data)
mac_mouse < subset(GSE146613,idents=c("Immunostimulatory_Mac","Immunosuppressive_Mac"))
my_comparisons <- list( c("Immunostimulatory_Mac","Immunosuppressive_Mac"))
p10 <- VlnPlot(mac, features = "Hpgds") +
    NoLegend()  +
    stat_compare_means(comparisons = my_comparisons)+ylim(-0.1,3.5)+
    scale_x_discrete(labels = c("Immunostimulatory_Mac" = "Immunostimulatory","Immunosuppressive_Mac" = "Immunosuppressive"))+
    theme(legend.position = "", 
          legend.title = element_text(face="bold", size=20),
          legend.text = element_text(size=20),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 15, face = "bold", angle = 0,hjust = 0.5,vjust = 1),        
          axis.title = element_text(face="bold", size=22),
          axis.title.x = element_text(color = "black", size = 15, face="bold"),
          plot.title = element_text(size = rel(2.5), face ="bold", hjust = 0.5,
                                    margin = margin(t = 10, b = 20, unit ="pt")))
ggsave(p10,filename = "HPGDS_M1vsM2.jpg",width = 9,height = 9)


###Figure S1A###
GC_immune <- GC_immune %>% subset(NEW_cluster_all_1 != "Melanocytes")
GC_immune <- GC_immune %>% subset(NEW_cluster_all_1 != "Doublets")
### FeaturePlot of HPGDS, all cells. Further aesthetical improvements were done over it but the base is the standard Seurat FeaturePlot
FeaturePlot(GC_immune, features = 'HPGDS', label = T)

# Proportion plot 
macrophages <- GC_immune %>% subset(NEW_cluster_all_1 %in% c( "M2 Macrophages",
                                                              'Macrophages_CXCL10',
                                                              'Macrophages_FOS',
                                                              'Macrophages_necrosis',
                                                              'Macrophages_PLA2G2D'))

macrophages@meta.data <- macrophages@meta.data %>%
  mutate(NEW_cluster_all_1 = recode(NEW_cluster_all_1 , "Macrophages_necrosis" = "Macrophages_SPP1",
                                                        "Macrophages_CXCL10" = "Macrophages_CXCL9",
                                                        "Macrophages_FOS" = "Macrophages_CCL3",
                                                        "M2 Macrophages"= "Macrophages_LYVE1"))

# take the metadata
TEST <- macrophages@meta.data
TEST$NEW_cluster_all_1<-gsub(x = TEST$NEW_cluster_all_1,  pattern = '_',replacement = '.')
TEST$NEW_cluster_all_1<-gsub(x = TEST$NEW_cluster_all_1,  pattern = ' ',replacement = '.')
levels(as.factor(TEST$NEW_cluster_all_1))
TEST$gc_number<-TEST$`GC number`
# make percentages
cell_num <- TEST %>%
  mutate(sample_id = as.factor(paste(NEW_cluster_all_1, `Response?`,`BT/OT`, orig.ident, sep="_"))) %>%
  mutate(gc_number = as.factor(gc_number)) %>%
  group_by(sample_id, gc_number, .drop=FALSE) %>%  dplyr::summarise(n=n()) %>%
  tidyr::separate(sample_id, c("NEW_cluster_all_1", "Response", "Timepoint", "orig.ident"),sep = '_')
cell_num

total_cells<- TEST %>%
  group_by(NEW_cluster_all_1) %>%
  dplyr::summarise(total = n())
total_cells

cell_percentage<- left_join(cell_num, total_cells) %>%
  mutate(percentage = n/total*100)
cell_percentage

cell_percentage <- subset(cell_percentage, Response != "NA")

cell_percentage <- cell_percentage %>%
  mutate(Response = recode(Response , "0" = "NR", "1" = "R"))

cell_percentage$var <- paste(cell_percentage$Timepoint,cell_percentage$Response, sep = "_" )

cell_percentage <- cell_percentage %>%
  mutate(gc_number = ifelse(orig.ident == "scrCMA040", "16_bis", gc_number))
cell_percentage <- cell_percentage %>%
  mutate(gc_number = ifelse(orig.ident == "scrCMA048", "16_bis", gc_number))

cell_percentage$patient_number<-cell_percentage$gc_number

p11 <- ggplot(cell_percentage, aes(NEW_cluster_all_1, 
                                y=percentage, 
                                fill=patient_number)) + 
  geom_bar(stat="identity", position = "stack")  + 
  theme_classic() + RotatedAxis() +  theme(text=element_text(size=12), 
                                           axis.text.x = element_text(angle = 45,hjust = 1),
                                           plot.title = element_text(hjust = 0.5))
ggsave(p11, file="percentage_macrophages.jpg", width = 8, height = 6)


###Figure S1D###
object.markers <- FindMarkers(OT, ident.1 = "NR", ident.2 = "R", group.by = "Response", logfc.threshold = 0, min.pct = 0, pseudocount.use=0.01)
cluster5.markers <- FindMarkers(mac, ident.1="Macrophages_LYVE1", ident.2 = c("Macrophages_PLA2G2D","Macrophages_CCL3","Macrophages_CXCL9","Macrophages_SPP1"), min.pct = 0.25)
marker <- rownames(cluster5.markers)
object <- object.markers[rownames(object.markers) %in% marker,]    ##(gene list)
cluster1.markers <- object %>%
     mutate(Difference = pct.1 - pct.2) %>% 
     rownames_to_column("gene")
 
cluster1.markers$threshold="ns"
cluster1.markers[which(cluster1.markers$avg_log2FC  >= 1.5 & cluster1.markers$p_val_adj <0.05),]$threshold="up";
cluster1.markers[which(cluster1.markers$avg_log2FC  <= -1.5 & cluster1.markers$p_val_adj < 0.05),]$threshold="down"

p14 <- ggplot(cluster1.markers, aes(x=Difference, y=avg_log2FC, color=threshold)) + 
     geom_point(size=2) + 
     scale_color_manual(name="differential expression",
                        values=c("grey","red")) + 
     geom_label_repel(data=subset(cluster1.markers, avg_log2FC >= 1.5 &  p_val_adj <= 0.05), 
                      aes(label=gene),  
                      color="black",
                      segment.colour = "black",
                      segment.size = 0.5,  
                      size=7)+
     geom_vline(xintercept = 0.0,linetype=4)+
     geom_hline(yintercept = 0,linetype=4)+
     theme_classic()+theme(axis.line = element_line(color = "black", size = 1),
                           axis.text.y = element_text(size = 12),
                           axis.text.x = element_text(size = 12,hjust = 1,vjust = 1))
ggsave(p14,filename = "volcano_m2_marker.jpg",width = 8,height = 6)


###Figure S1E###
mac@meta.data <- tidyr::unite(mac@meta.data,"RT", "Response", "BT/OT",sep="_",remove=FALSE)
df <- VlnPlot(mac, features = "HPGDS", group.by="RT",split.by = "NEW_cluster_all_1")$data
df$split<- factor(df$split, levels=c("Macrophages_LYVE1","Macrophages_PLA2G2D","Macrophages_SPP1","Macrophages_CCL3","Macrophages_CXCL9"))

stat.test <- compare_means(
     HPGDS ~ ident, data = df, group.by = "split",
     method = "wilcox.test", ref.group = "NR_OT"
 )
stat.test <- stat.test[-c(1,2),]
stat.test <- stat.test[-c(2:13),]
stat.test <- stat.test %>%
     mutate(y.position = c(2))

p15 <- ggplot(data = df,aes(x=ident,y = HPGDS)) +
    geom_violin(aes(fill = ident), trim=TRUE,scale = "width")+
    scale_y_continuous(limits =  c(-0.01,3),position="left",labels = function(x)
        c(rep(x = "", times = length(x)-2), x[length(x) - 1], ""))+
    facet_grid(rows = vars(split), scales = 'free_x')+ scale_x_discrete(limits = c("NR_BT","R_BT","NR_OT","R_OT"))+
    theme(legend.position = "none", panel.spacing = unit(0, "lines"),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = NA, color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.text.y.right = element_text(angle = 0, size = 12, hjust = 0,face = "bold", vjust = 0.5, color = "black"),
          strip.text.y.left = element_text(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 45,size = 12, hjust = 0.5,face = "bold", vjust = 0.5, color = "black")
    )+scale_fill_manual(values=c("#F8766D","#F8766D","#00B8E7","#00B8E7"))
p16 <- p15+stat_pvalue_manual(stat.test, label = "p.adj", size=3)
p17 <- p16 + theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

ggsave(p17 ,file="vln1.jpg",height = 3, width = 5)


###Figure S1F###
GC_all_immune <- readRDS('GC_all_immune_final_version.rds')

table(GC_all_immune$high_resolution_clusters)
macrophages <- GC_all_immune %>% subset(high_resolution_clusters %in% c("Macrophages_CXCL9",
                                                                        "Macrophages_LYVE1",
                                                                        "Macrophages_CCL3",
                                                                        "Macrophages_SPP1",
                                                                        "Macrophages_PLA2G2D")
)
levels(as.factor(macrophages$high_resolution_clusters))
macrophages$HPGDS_positive <- ifelse(macrophages@assays$SCT@data["HPGDS",] > 0,
                                     "HPGDS_positive",
                                     "HPGDS_negative")
#
TEST <- macrophages@meta.data
TEST$high_resolution_clusters <- gsub(x = TEST$high_resolution_clusters,
                                      pattern = '_',
                                      replacement = '.')
levels(as.factor(TEST$high_resolution_clusters))
#
# make percentages
cell_num <- TEST %>%
  mutate(sample_id = as.factor(paste(sample_ID,
                                     Response,
                                     Timepoint, 
                                     patient_ID,
                                     sep="_"))) %>%
  mutate(HPGDS_positive = as.factor(HPGDS_positive)) %>%
  group_by(sample_id, HPGDS_positive, .drop=FALSE) %>%  dplyr::summarise(n=n()) %>%
  tidyr::separate(sample_id, c("sample_ID", 
                               "Response", 
                               "Timepoint", 
                               "patient_ID"))
cell_num
total_cells<- TEST %>%
  group_by(sample_ID) %>%
  dplyr::summarise(total = n())
total_cells

cell_percentage<- left_join(cell_num, total_cells) %>%
  mutate(percentage = n/total*100)
cell_percentage
cell_percentage$var <- paste(cell_percentage$Timepoint,cell_percentage$Response, sep = "_")
cell_percentage$var <- ordered(cell_percentage$var, levels = c("BT_NR", "OT_NR", "BT_R",'OT_R'))

# keep only positive for cell_percentage 
cell_percentage <- cell_percentage[cell_percentage$HPGDS_positive=='HPGDS_positive',]
cell_percentage

# Rename sample
cell_percentage <- cell_percentage %>%
  mutate(patient_ID = ifelse(sample_ID == "scrCMA040", "16_bis", patient_ID))
cell_percentage <- cell_percentage %>%
  mutate(patient_ID = ifelse(sample_ID == "scrCMA048", "16_bis", patient_ID))

comparison1 <- list(c("BT_NR", "OT_NR"), c("BT_R",'OT_R'))

p18 <- ggboxplot(cell_percentage, x = "var", y = "percentage",
          fill = "Timepoint",
          shape = "Response",
          palette = c("#F8766D","#00B8E7"),
          add = "jitter",
          order = c("BT_NR", "OT_NR", "BT_R",'OT_R')) +
          #
          geom_signif(comparisons = comparison1, 
                      step_increase = 0.1, 
                      test = 'wilcox.test', 
                      test.args = list(paired=T)) +
          #
          geom_line(aes(group = patient_ID), alpha = 0.6, colour = "black") +
          #  
          facet_wrap(~HPGDS_positive, scales = "free_y") +
          theme_classic() + 
          RotatedAxis() +  
          theme(text=element_text(size=12),
                axis.text.x = element_text(angle = 90,hjust = 1),
                plot.title = element_text(hjust = 0.5)
                )
ggsave(p18, file="paired analysis.jpg", height = 4, width = 8)

                       
###Figure S1G###
Idents(GSE120575)=GSE120575$celltype
mac_cluster <- subset(GSE120575, idents = c("Monocyte", "Macrophage"))
view(mac_cluster@meta.data)
mac_cluster@meta.data <- separate(mac_cluster@meta.data, col = "characteristics..patinet.ID..Pre.baseline..Post..on.treatment.", into = c("treat","postID"), sep = '_')
mac_cluster@meta.data <- tidyr::unite(mac_cluster@meta.data,"Cell",treat,characteristics..response,sep="_",remove=FALSE)
my_comp <- list(c("Post_Non-responder", "Post_Responder"),c("Pre_Non-responder", "Pre_Responder"))
Idents(mac_cluster) <- mac_cluster$characteristics..therapy
mac_cluster1 <- subset(mac_cluster, idents = c("anti-PD1", "anti-CTLA4+PD1"))
vln_df = data.frame(HPGDS = mac_cluster1[["RNA"]]@data["HPGDS",], cluster = mac_cluster1$Cell)
p19 <- ggplot(vln_df, aes(x = cluster, y = HPGDS)) + geom_violin(aes(fill = cluster), trim=TRUE, scale = "width") + 
    stat_compare_means(comparisons = my_comp, method = "wilcox.test", label = "p.adj", size=3)+ylim(-0,15)+
    theme_classic() +
    ggtitle("")+
    theme(legend.position = "", 
          legend.title = element_text(face="bold", size=20),
          legend.text = element_text(size=20),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 8, face = "bold",color = "black", angle = 45,hjust = 1,vjust = 1),        
          axis.title = element_text(face="bold", size=15),
          axis.title.x = element_text(color = "black", size = 15, face="bold"),
          plot.title = element_text(size = rel(2.5), face ="bold", hjust = 0.5,
                                    margin = margin(t = 10, b = 20, unit ="pt")))+
    xlab("")+scale_x_discrete(limits = c("Pre_Non-responder","Pre_Responder","Post_Non-responder","Post_Responder"))+
    scale_fill_manual(values=c("#F8766D","#00B8E7","#F8766D","#00B8E7"))

ggsave(p19, file="HPGDS.jpg", height = 4, width = 5)
