##### Bulk RNA-seq of Hpgds and WT KO TAMs - R Source code #####

#Set seed for reproducibility
set.seed(1234)

#Load necessary packages
library("DESeq2")
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)
library(ggtext)
library(ggrepel)
library(xlsx)
library(stringr)
library(ggforce)
library(clusterProfiler)
library(enrichplot)
library(plyr)
library("org.Mm.eg.db", character.only = TRUE)
require(DOSE)
library(GSVA)
library(msigdbr)
library(fgsea)
library(readxl)

### Set working directory in results folder 
setwd("/results")
getwd()

### Load countmatrix
countData <- as.matrix(as.data.frame(utils::read.csv("/data/genecounts_s2.counts.counts", sep="\t", row.names="Geneid")))

### Prepare metadata
    #Simplify colnames; remove everything after first underscore
    colnames(countData) <- gsub("_.*","", colnames(countData))
    #Define conditions and compile all metadata into dataframe
    revalue(colnames(countData), c(
        'RT1' = "WT", 'RT2' = "WT", 'RT4' = "WT",
        'RT6' = "Hpgds KO", 'RT7' = "Hpgds KO", 'RT8' = "Hpgds KO")) -> condition
    condition <- factor(condition, levels = c("WT", "Hpgds KO"))
    colData <- data.frame(sampleName = colnames(countData), condition = condition)

### Create DESEQ2 object and perform normalization
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = colData, design =~ condition)
    colnames(dds) <- colnames(countData) #ensure the column names are correct.
    dds <- DESeq2::DESeq(dds, full =~ condition, fitType = "local")

### Compute DEGs with shrunken log2 fold changes with apeglm
    # If needed, confirm which coef to select: 
    DESeq2::resultsNames(dds)
    #Extract DEG information
    res_shrink <- DESeq2::lfcShrink(dds, coef = 2, type = "apeglm")
    #Remove genes with missing adjusted p-values
    res_shrinkFix <- res_shrink[!is.na(res_shrink$padj),]
    #Convert to dataframe
    res_shrinkDF <- as.data.frame(res_shrinkFix)
    #Save output as xlsx files
    xlsx::write.xlsx(res_shrinkDF, "DEG_Result_Wald_Test.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)

### Visualize DEG analyses as a Volcanoplot and save as pdf file
    #Add column to determine dot color
    up_inKO <- (res_shrinkDF$log2FoldChange > 0 & res_shrinkDF$padj < 0.05)
    up_inWT <- (res_shrinkDF$log2FoldChange < 0 & res_shrinkDF$padj < 0.05)
    res_shrinkDF$DE <- ifelse(up_inWT,'up in WT',ifelse(up_inKO,'up in Hpgds KO','not DE'))
    #Plot actual graph
    Graph <- ggplot2::ggplot(res_shrinkDF,aes(log2FoldChange,-log10(pvalue))) + 
        geom_point(aes(color=DE)) +
        scale_color_manual(name="differential expression",
                     values=c("black", "#FD8008","#999999"),
                     labels=str_wrap(c("Not significant", "Up in Hpgds KO TAMs", "Down in Hpgds KO TAMs"), width = 45)) + 
        ylab("-log10 p-value") + 
        xlab("log2 fold change") + 
        theme_classic() +
        theme(axis.text.x = element_text(color="black", size = 12), axis.text.y = element_text(color="black", size = 12),
              axis.title.x = element_text(color="black", size = 15), axis.title.y = element_text(color="black", size = 15),
              aspect.ratio = 1) +
        ggrepel::geom_text_repel(force=5, col= "red", aes(label=ifelse(rownames(res_shrinkDF) == "Hpgds", as.character(rownames(res_shrinkDF)),'')), max.overlaps = Inf) + #Label Hpgds specifically
        ggrepel::geom_text_repel(force=50, col= "black", aes(label=ifelse(rownames(res_shrinkDF) %in% c("Rims3", "Tmc3", "2010016I18Rik", "Myom1", "Scarf2", "Negr1","Adam22",
                                                                                                       "Ly6a2", "Ltc4s", "Ly6c2", "Ctnnd2", "Xlr3b", "Ccl2"), as.character(rownames(res_shrinkDF)),'')), max.overlaps = Inf) #Label other genes specifically
    # Save as PDF
    pdf("VolcanoPlot_DEG_Hpgds_Labelled.pdf", useDingbats = FALSE)
        plot(Graph)
    dev.off()   

### Visualize the expression of Hpgds in WT and Hpgds KO TAMs.
    # Getting the count table in ggplot-friendly format
    geneCounts = DESeq2::plotCounts(dds, gene="Hpgds", intgroup=c("condition"), returnData=TRUE)
    #Plot Graph using ggplot2
    Graph <- ggplot2::ggplot(geneCounts, aes(x=condition, y=count, fill = condition, col = condition)) +
         geom_bar(position = "dodge", stat = "summary") +
         scale_fill_manual(values = c("WT" = "#D4D4D4", "Hpgds KO" = "#FFE0C0")) + 
         scale_colour_manual(values = c("WT" = "black", "Hpgds KO" = "#FF5700")) + 
         geom_point(aes(col = condition), size = 3, shape = 15) +
         theme_classic()+
         theme(axis.text.x = element_text(color="black", size = 15), 
               axis.text.y = element_text(color="black", size = 15),
               legend.position = "none",
               aspect.ratio = 1.4,
               plot.title = element_markdown(size = 15, hjust = 0.5),
               axis.title = element_text(size=15)) +
         ggtitle('*Hpgds* expression in murine TAMs')+
         ylab('Normalized counts')+
         xlab('') +
         scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
    # Save graph as pdf file
    pdf("Barplot_Hpgds.pdf", useDingbats = FALSE)
        plot(Graph)
    dev.off()   

### Visualize the expression of Hgf in WT and Hpgds KO TAMs.
    # Getting the count table in ggplot-friendly format
    geneCounts = DESeq2::plotCounts(dds, gene="Hgf", intgroup=c("condition"), returnData=TRUE)
    #Plot Graph using ggplot2
    Graph <- ggplot2::ggplot(geneCounts, aes(x=condition, y=count, fill = condition, col = condition)) +
         geom_bar(position = "dodge", stat = "summary") +
         scale_fill_manual(values = c("WT" = "#D4D4D4", "Hpgds KO" = "#FFE0C0")) + 
         scale_colour_manual(values = c("WT" = "black", "Hpgds KO" = "#FF5700")) + 
         geom_point(aes(col = condition), size = 3, shape = 15) +
         theme_classic()+
         theme(axis.text.x = element_text(color="black", size = 15), 
               axis.text.y = element_text(color="black", size = 15),
               legend.position = "none",
               aspect.ratio = 1.4,
               plot.title = element_markdown(size = 15, hjust = 0.5),
               axis.title = element_text(size=15)) +
         ggtitle('*Hgf* expression in murine TAMs')+
         ylab('Normalized counts')+
         xlab('') +
         scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
    # Save graph as pdf file
    pdf("Barplot_Hgf.pdf", useDingbats = FALSE)
        plot(Graph)
    dev.off() 

### Perform GSEA
    #Prepare ranked genelist, sorted by shrunken log2 Fold Change, for GSEA 
    original_gene_list <- res_shrinkFix$log2FoldChange
    names(original_gene_list) <- rownames(res_shrinkFix)
    gene_list <- stats::na.omit(original_gene_list) # omit any NA values 
    gene_list = sort(gene_list, decreasing = TRUE) # sort the list in decreasing order
    # Perform GSEA
    ### Biological Process ###
    gse_BP <- clusterProfiler::gseGO(geneList = gene_list, 
              ont ="BP", 
              keyType = "SYMBOL", #We use gene ID here
              pvalueCutoff = 1, #We want to retrieve whole output 
              verbose = TRUE, 
              OrgDb = "org.Mm.eg.db", 
              pAdjustMethod = "fdr", 
              eps = 0,
              seed = TRUE) #Set seed for reproducibility
    ### Molecular Function ###
    gse_MF <- clusterProfiler::gseGO(geneList = gene_list, 
              ont ="MF", 
              keyType = "SYMBOL", #We use gene ID here
              pvalueCutoff = 1, #We want to retrieve whole output 
              verbose = TRUE, 
              OrgDb = "org.Mm.eg.db", 
              pAdjustMethod = "fdr", 
              eps = 0,
              seed = TRUE) #Set seed for reproducibility
    ### Cellular Component ###
    gse_CC <- clusterProfiler::gseGO(geneList = gene_list, 
              ont ="CC",
              keyType = "SYMBOL", #We use gene ID here
              pvalueCutoff = 1, #We want to retrieve whole output 
              verbose = TRUE, 
              OrgDb = "org.Mm.eg.db", 
              pAdjustMethod = "fdr", 
              eps = 0,
              seed = TRUE) #Set seed for reproducibility
    
    #Save output as xlsx files
    ### Biological Process ###
    results <- as.data.frame(gse_BP)
    xlsx::write.xlsx(results, "GSEA_BiologicalProcess_completeresults.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)
    ### Molecular Function ###
    results <- as.data.frame(gse_MF)
    xlsx::write.xlsx(results, "GSEA_MolecularFunction_completeresults.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)
    ### Cellular Component ###
    results <- as.data.frame(gse_CC)
    xlsx::write.xlsx(results, "GSEA_CellularComponent_completeresults.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)


        #Generate dotplot with significant (adjusted p-value < 0.05) GSEA results and save as PDF
    ### Biological Process ###
    Graph <- dotplot(gse_BP %>% filter(p.adjust < 0.05), 
                     showCategory=450, 
                     split=".sign",       
                     label_format = 320) + 
                     facet_grid(.~.sign) +
             ggtitle("GSEA - Biological Process") +
             theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))
    pdf("GSEA_Biological_Process.pdf", useDingbats = FALSE, width = 10, height = 10)
    plot(Graph)
    dev.off() 
    ### Molecular Function ###
        Graph <- dotplot(gse_MF %>% filter(p.adjust < 0.05), 
                     showCategory=400, 
                     split=".sign",       
                     label_format = 320) + 
                     facet_grid(.~.sign) +
             ggtitle("GSEA - Molecular Function") +
             theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))
    pdf("GSEA_Molecular_Function.pdf", useDingbats = FALSE, width =10, height = 8)
    plot(Graph)
    dev.off() 
    ### Cellular Component ###
        Graph <- dotplot(gse_CC %>% filter(p.adjust < 0.05), 
                     showCategory=450, 
                     split=".sign",       
                     label_format = 320) + 
                     facet_grid(.~.sign) +
             ggtitle("GSEA - Cellular Component") +
             theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))
    pdf("GSEA_Cellular_Component.pdf", useDingbats = FALSE, width = 10, height = 8)
    plot(Graph)
    dev.off() 

### Perform fGSEA with custom gene signatures

# DEG induced/repressed in murine macrophages by different cytokines (Available as Table S5 of the Cui,Hacogen et al. Nature 2023 paper (https://www.nature.com/articles/s41586-023-06816-9#MOESM3))
    Cytokine_DEG <- readxl::read_xlsx("/data/41586_2023_6816_MOESM5_ESM.xlsx", sheet = "Macrophage")
    #Extract DEG from macrophages induced by different cytokines
    Cytokine_DEG_Macrophage <- Cytokine_DEG %>% filter(Celltype_Str == "Macrophage") %>% split(Cytokine_DEG, f = Cytokine_DEG$Cytokine_Str)
    #Prepare genesets for GSEA
    custom_gsea_sets <- list(
    "Response to Cardiotrophin1" = unique(Cytokine_DEG_Macrophage$`Cardiotrophin-1`$Gene),
    "Response to Decorin" = unique(Cytokine_DEG_Macrophage$Decorin$Gene),
    "Response to EGF" = unique(Cytokine_DEG_Macrophage$EGF$Gene),
    "Response to GM-CSF" = unique(Cytokine_DEG_Macrophage$`GM-CSF`$Gene),
    "Response to IFNa1" = unique(Cytokine_DEG_Macrophage$IFNa1$Gene),
    "Response to IFNb" = unique(Cytokine_DEG_Macrophage$IFNb$Gene),
    "Response to IFNe" = unique(Cytokine_DEG_Macrophage$IFNe$Gene),
    "Response to IFNg" = unique(Cytokine_DEG_Macrophage$IFNg$Gene),
    "Response to IFNk" = unique(Cytokine_DEG_Macrophage$IFNk$Gene),
    "Response to IL10" = unique(Cytokine_DEG_Macrophage$IL10$Gene),
    "Response to IL11" = unique(Cytokine_DEG_Macrophage$IL11 $Gene),
    "Response to IL12" = unique(Cytokine_DEG_Macrophage$IL12$Gene),
    "Response to IL13" = unique(Cytokine_DEG_Macrophage$IL13$Gene),
    "Response to IL15" = unique(Cytokine_DEG_Macrophage$IL15$Gene),
    "Response to IL17F" = unique(Cytokine_DEG_Macrophage$IL17F$Gene),
    "Response to IL18" = unique(Cytokine_DEG_Macrophage$IL18$Gene),
    "Response to IL1a" = unique(Cytokine_DEG_Macrophage$IL1a$Gene),
    "Response to IL1b" = unique(Cytokine_DEG_Macrophage$IL1b$Gene),
    "Response to IL2" = unique(Cytokine_DEG_Macrophage$IL2$Gene),
    "Response to IL21" = unique(Cytokine_DEG_Macrophage$IL21$Gene),
    "Response to IL24" = unique(Cytokine_DEG_Macrophage$IL24$Gene),
    "Response to IL27" = unique(Cytokine_DEG_Macrophage$IL27$Gene),
    "Response to IL3" = unique(Cytokine_DEG_Macrophage$IL3$Gene),
    "Response to IL31" = unique(Cytokine_DEG_Macrophage$IL31$Gene),
    "Response to IL33" = unique(Cytokine_DEG_Macrophage$IL33$Gene),
    "Response to IL36a" = unique(Cytokine_DEG_Macrophage$IL36a$Gene),
    "Response to IL36RA" = unique(Cytokine_DEG_Macrophage$IL36RA$Gene),
    "Response to IL4" = unique(Cytokine_DEG_Macrophage$IL4$Gene),
    "Response to IL5" = unique(Cytokine_DEG_Macrophage$IL5$Gene),
    "Response to IL7" = unique(Cytokine_DEG_Macrophage$IL7$Gene),
    "Response to Leptin" = unique(Cytokine_DEG_Macrophage$Leptin$Gene),
    "Response to LIF" = unique(Cytokine_DEG_Macrophage$LIF$Gene),
    "Response to M-CSF" = unique(Cytokine_DEG_Macrophage$`M-CSF`$Gene),
    "Response to OSM" = unique(Cytokine_DEG_Macrophage$OSM$Gene),
    "Response to Persephin" = unique(Cytokine_DEG_Macrophage$Persephin$Gene),
    "Response to Prolactin" = unique(Cytokine_DEG_Macrophage$Prolactin$Gene),
    "Response to RANKL" = unique(Cytokine_DEG_Macrophage$RANKL$Gene),
    "Response to SCF" = unique(Cytokine_DEG_Macrophage$SCF$Gene),
    "Response to TGF-beta-1" = unique(Cytokine_DEG_Macrophage$`TGF-beta-1`$Gene),
    "Response to TNFa" = unique(Cytokine_DEG_Macrophage$TNFa$Gene),
    "Response to TPO" = unique(Cytokine_DEG_Macrophage$TPO$Gene))

    #Prepare ranked genelist as input
    res_shrinkDF$feature <- rownames(res_shrinkDF)
    rankedDEG <- res_shrinkDF %>%
                              dplyr::arrange(desc(log2FoldChange)) %>%
                              dplyr::select(feature, log2FoldChange)
    ranks <- tibble::deframe(rankedDEG)
    #Perform GSEA itself
    fgseaRes<- fgsea::fgsea(custom_gsea_sets, stats = ranks)
    #Tidy output
    fgseaResTidy <- fgseaRes %>%
        tibble::as_tibble() %>%
        dplyr::arrange(desc(NES))
    #Convert column from list to string to allow saving of output as xlsx file
    fgseaResTidy$leadingEdge  <-  sapply(fgseaResTidy$leadingEdge, toString)
    #Save output as xlsx file
    xlsx::write.xlsx(fgseaResTidy, "fGSEA_CytokineDictionaryMacrophages.xlsx", sheetName = "Sheet1", append = FALSE, row.names = TRUE)
    #Plot a barplot for with the normalized Enrichment score and save as pdf
    Graph <- ggplot2::ggplot(fgseaResTidy, aes(NES,reorder(pathway, NES))) +
                geom_col(aes(fill= NES < 0)) +
                geom_text(aes(label = formatC(padj, format = "e")), position = position_stack(vjust = 0.5), size= 4)+
                labs(y="", x="Normalized Enrichment Score (fGSEA)", title="Macrophage Response to Cytokines")  +
                geom_vline(xintercept = 0, linetype = "solid", color = "black") +
                coord_cartesian(xlim=c(-2,2)) +
                theme_classic() +
                theme(aspect.ratio=3, legend.position = "none",
                      axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black", size = 12),
                      plot.title = element_text(hjust = 0.5, face = "bold"),
                      panel.grid.major.y = element_line(), panel.grid.minor.y = element_line()) +
                scale_fill_manual( values = c("#FD8008","#999999")) #Custom colors barplot
    
    pdf("fGSEA_barplot_CytokineDictionaryMacrophages.pdf", height = 10, width = 10, useDingbats = F)
        plot(Graph)
    dev.off()

# DEG linked to different macrophage polarization states (Available as Table S7 of the Cui,Hacogen et al. Nature 2023 paper (https://www.nature.com/articles/s41586-023-06816-9#MOESM3))
    MacrophagePolarization_DEG <- readxl::read_xlsx("/data/41586_2023_6816_MOESM9_ESM.xlsx", sheet = "Macrophage")
    #Extract DEG from different macrophage polarization states
    MacrophagePolarization_DEG <- MacrophagePolarization_DEG %>% split(MacrophagePolarization_DEG, f = MacrophagePolarization_DEG$Polarization)
    #Prepare genesets for GSEA
    custom_gsea_sets <- list(
    "Macrophage polarization state a" = unique(MacrophagePolarization_DEG$`Mac-a`$Gene),
    "Macrophage polarization state b" = unique(MacrophagePolarization_DEG$`Mac-b`$Gene),
    "Macrophage polarization state c" = unique(MacrophagePolarization_DEG$`Mac-c`$Gene),
    "Macrophage polarization state d" = unique(MacrophagePolarization_DEG$`Mac-d`$Gene),
    "Macrophage polarization state e" = unique(MacrophagePolarization_DEG$`Mac-e`$Gene))
    #Perform GSEA itself
    fgseaRes<- fgsea::fgsea(custom_gsea_sets, stats = ranks)
    #Tidy output
    fgseaResTidy <- fgseaRes %>%
        tibble::as_tibble() %>%
        dplyr::arrange(desc(NES))
    #Convert column from list to string to allow saving of output as xlsx file
    fgseaResTidy$leadingEdge  <-  sapply(fgseaResTidy$leadingEdge, toString)
    #Save output as xlsx file
    xlsx::write.xlsx(fgseaResTidy, "fGSEA_MacrophagePolarizationState.xlsx", sheetName = "Sheet1", append = FALSE, row.names = TRUE)
    #Plot a barplot for with the normalized Enrichment score and save as pdf
    Graph <- ggplot2::ggplot(fgseaResTidy, aes(NES,reorder(pathway, NES))) +
                geom_col(aes(fill= NES < 0)) +
                geom_text(aes(label = formatC(padj, format = "e")), position = position_stack(vjust = 0.5), size= 4)+
                labs(y="", x="Normalized Enrichment Score (fGSEA)", title="Macrophage polarization state")  +
                geom_vline(xintercept = 0, linetype = "solid", color = "black") +
                coord_cartesian(xlim=c(-2,2)) +
                theme_classic() +
                theme(aspect.ratio=2, legend.position = "none",
                      axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black", size = 12),
                      plot.title = element_text(hjust = 0.5, face = "bold"),
                      panel.grid.major.y = element_line(), panel.grid.minor.y = element_line()) +
                scale_fill_manual( values = c("#FD8008","#999999")) #Custom colors barplot
    
    pdf("fGSEA_barplot_MacrophagePolarizationState.pdf", height = 5, width = 10, useDingbats = F)
        plot(Graph)
    dev.off()

# M1 M2 macrophage genesets from this paper
    custom_gsea_sets <- list(
        "Mrc1_immunosuppressive_macrophages" = unique(c("C1qa","Mrc1","Sepp1","Ctsb","Retnla","Cd200r1")))
    #Perform GSEA itself
    fgseaRes<- fgsea::fgsea(custom_gsea_sets, stats = ranks)
    #Tidy output
    fgseaResTidy <- fgseaRes %>%
        tibble::as_tibble() %>%
        dplyr::arrange(desc(NES))
    #Convert column from list to string to allow saving of output as xlsx file
    fgseaResTidy$leadingEdge  <-  sapply(fgseaResTidy$leadingEdge, toString)
    #Save output as xlsx file
    xlsx::write.xlsx(fgseaResTidy, "fGSEA_M1_M2_Macrophage.xlsx", sheetName = "Sheet1", append = FALSE, row.names = TRUE)
    #Plot a barplot for with the normalized Enrichment score and save as pdf
    Graph <- ggplot2::ggplot(fgseaResTidy, aes(NES,reorder(pathway, NES))) +
                geom_col(aes(fill= NES < 0)) +
                geom_text(aes(label = formatC(padj, format = "e")), position = position_stack(vjust = 0.5), size= 4)+
                labs(y="", x="Normalized Enrichment Score (fGSEA)", title="")  +
                geom_vline(xintercept = 0, linetype = "solid", color = "black") +
                coord_cartesian(xlim=c(-2,2)) +
                theme_classic() +
                theme(aspect.ratio=2, legend.position = "none",
                      axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black", size = 12),
                      plot.title = element_text(hjust = 0.5, face = "bold"),
                      panel.grid.major.y = element_line(), panel.grid.minor.y = element_line()) +
                scale_fill_manual( values = c("#FD8008","#999999")) #Custom colors barplot
    
    pdf("fGSEA_barplot_M1_M2_Macrophage.pdf", height = 5, width = 10, useDingbats = F)
        plot(Graph)
    dev.off()

# Angiogenesis signature
    custom_gsea_sets <- list(
        "Angiogenesis" = unique(c("Fgf2","Pdgfc","Egf","Hgf","Pdgfa","Tgfa","Angpt2","Vegfc","Angpt1","Pdgfb","Vegfa","Vegfb", "Igf1", "Pigf")))
    #Perform GSEA itself
    fgseaRes<- fgsea::fgsea(custom_gsea_sets, stats = ranks)
    #Tidy output
    fgseaResTidy <- fgseaRes %>%
        tibble::as_tibble() %>%
        dplyr::arrange(desc(NES))
    #Convert column from list to string to allow saving of output as xlsx file
    fgseaResTidy$leadingEdge  <-  sapply(fgseaResTidy$leadingEdge, toString)
    #Save output as xlsx file
    xlsx::write.xlsx(fgseaResTidy, "fGSEA_Angiogenesis.xlsx", sheetName = "Sheet1", append = FALSE, row.names = TRUE)
    #Plot a barplot for with the normalized Enrichment score and save as pdf
    Graph <- ggplot2::ggplot(fgseaResTidy, aes(NES,reorder(pathway, NES))) +
                geom_col(aes(fill= NES < 0)) +
                geom_text(aes(label = formatC(padj, format = "e")), position = position_stack(vjust = 0.5), size= 4)+
                labs(y="", x="Normalized Enrichment Score (fGSEA)", title="")  +
                geom_vline(xintercept = 0, linetype = "solid", color = "black") +
                coord_cartesian(xlim=c(-2,2)) +
                theme_classic() +
                theme(aspect.ratio=2, legend.position = "none",
                      axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black", size = 12),
                      plot.title = element_text(hjust = 0.5, face = "bold"),
                      panel.grid.major.y = element_line(), panel.grid.minor.y = element_line()) +
                scale_fill_manual( values = c("#FD8008","#999999")) #Custom colors barplot
    
    pdf("fGSEA_barplot_Angiogenesis.pdf", height = 5, width = 10, useDingbats = F)
        plot(Graph)
    dev.off()

# Show session information
    sessionInfo()

