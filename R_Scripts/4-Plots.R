library(Seurat)
library(Matrix)
library(dplyr)
library(dittoSeq)
library(scater)
library(patchwork)
library(readxl)
library(scuttle)
library(glmGamPoi)
library(sctransform)

setwd("E:/ASAP/scRNASeq/2D_DAN/all_WT_analysis/iSCORE-PD/")

all <- readRDS("all_WIBR3_DAN.rds")
foundin <- readRDS("../FOUNDIN-PD/iPSCsDopaALL_integratedAfterBroadCellType.RDS")
foundin <- RunUMAP(foundin, dims = 1:30, return.model=TRUE)

pseudoBulk_spearman <- function(obj_1, obj_2, labels_1, labels_2, assay, n_var_features){
    Idents(obj_1) <- "seurat_clusters"
    obj_1@active.assay <- assay
    obj_1 <- FindVariableFeatures(obj_1, assay = assay, nfeatures = n_var_features)
    Idents(obj_2) <- "CellType"
    obj_2@active.assay <- assay
    obj_2 <- FindVariableFeatures(obj_2, assay = assay, nfeatures = n_var_features)
    df <- data.frame(matrix(ncol = length(unique(obj_1$seurat_clusters)), nrow = length(unique(obj_2$CellType))))
    colnames(df) <- unique(obj_1$seurat_clusters)
    rownames(df) <- unique(obj_2$CellType)
    for (celltype in unique(obj_1$seurat_clusters)) {
        x <- (subset(obj_1, seurat_clusters == celltype))@assays[[assay]]@counts
        x_pseudobulk <- data.frame(Gene = rownames(x), MeanValue = rowMeans(x))
        x_pseudobulk <- x_pseudobulk[obj_1@assays[[assay]]@var.features,]
        for (celltype2 in unique(obj_2$CellType)) {
            y <- (subset(obj_2, CellType == celltype2))@assays[[assay]]@counts
            y_pseudobulk <- data.frame(Gene = rownames(y), MeanValue = rowMeans(y))
            y_pseudobulk <- y_pseudobulk[obj_2@assays[[assay]]@var.features,]
            gene_ids <- intersect(rownames(x_pseudobulk), rownames(y_pseudobulk))
            x_subset <- x_pseudobulk[gene_ids, , drop = FALSE]
            y_subset <- y_pseudobulk[gene_ids, , drop = FALSE]
            z <- cbind(x_subset, y_subset)
            spear <- cor(x = z[2], y = z[4], method = "spearman")
            df[celltype2,celltype] <- spear
        }
    }
    print(df)
    return(df)
}
markers <-function(x, phrase, var, min_pct, lfc_cutoff, test) {
    if (is.null(var)) {
        markers <- FindAllMarkers(x, only.pos = FALSE, min.pct = min_pct, logfc.threshold = lfc_cutoff, test.use = test) 
        markers %>%
            group_by(cluster) %>%
            top_n(n = 10, wt = avg_log2FC) -> top10
        p <- DoHeatmap(x, features = top10$gene) + NoLegend()
        ggsave(filename = paste0(phrase, "_MarkerHeatMap.png"), plot = p, device = "png", width = 30, height = 25)
        write.csv(markers, file = paste0(phrase, "_markers.csv"))
        write.csv(top10, file = paste0(phrase, "_top10markers.csv"))
    } else {
        x <- SetIdent(x, value = var)
        markers <- FindAllMarkers(x, only.pos = FALSE, min.pct = min_pct, logfc.threshold = lfc_cutoff, test.use = test) 
        markers %>%
            group_by(cluster) %>%
            top_n(n = 10, wt = avg_log2FC) -> top10
        p <- DoHeatmap(x, features = top10$gene) + NoLegend()
        ggsave(filename = paste0(phrase, "_MarkerHeatMap.png"), plot = p, device = "png", width = 30, height = 25)
        write.csv(markers, file = paste0(phrase, "_markers_by_", var, ".csv"))
        write.csv(top10, file = paste0(phrase, "_top10markers_by_", var, ".csv")) 
    }
}

####
#### Fig. 2C - UMAP with coarse cluster labels.
####
dittoDimPlot(all, "seurat_clusters_coarse", do.label = T, labels.size = 12, labels.repel = T)


####
#### Fig. 2D  - DotPlot of key marker genes
####
all@active.ident <- factor(all@active.ident, levels=c("0","1","4","6","3","5","2"))
goi <- c("TH","NR4A2","KCNJ6","LMX1B","LMX1A","STMN2",
           "MAP2","GAP43","SNAP25","SNCA","SYN1","ROBO1",
           "FOXA2","SLC18A1","SLC18A2","TPH1","MKI67",
           "HES1","SLIT2","FOXP2","CORIN","SOX6","CALB1")
DotPlot(all, features = goi, assay = 'SCT', dot.scale = 8, scale.by = "size") + scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "Oranges")) + RotatedAxis()


####
#### Fig. 2E  - Heatmap of Spearmann correlations of pseudobulked WIBR3 DAN clusters vs pseudobulked FOUNDIN-PD clusters.
####
mtx <- pseudoBulk_spearman(all, foundin, "seurat_clusters_coarse", "CellType", "RNA", n_var_features = 3000)
saveRDS(mtx, "all_v_FOUNDIN_FULL_Spearman_matrix_RNA3000.rds")
mtx <- readRDS("all_v_FOUNDIN_Spearman_matrix_RNA3000.rds")
pheatmap(mtx,
    main = "Heatmap of Spearman correlations \n WIBR3 vs. FOUNDIN-PD DA Neurons\n (RNA Assay, pseudobulked clusters) \n Input = 3000 variable genes", 
    fontsize = 12, 
    clustering_method = "ward.D2", 
    cutree_rows = 4, 
    cutree_cols = 3, 
    treeheight_row = 30, 
    treeheight_col = 30)


####
#### Fig. S3A  - UMAPs of WIBR3 DAN labeled with coarse cluster labels, and split by subclone identity.
####
dittoDimPlot(all, "seurat_clusters_coarse", split.by = "subclone_ID", do.label = T, labels.size = 12, labels.repel = T, size=1.8)


####
#### Fig. S3B  - Heatmap of top 10 cluster-defining markers in WIBR3 DAN dataset.
####
markers(all, "all_WIBR3_DAN_coarse_clusters_MAST", "seurat_clusters_coarse", 0.3, 0.5, "MAST")
m <- read.csv("all_WIBR3_DAN_coarse_clusters_MAST_markers_by_seurat_clusters_coarse.csv")
m %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(all, features = top10$gene, group.colors = dittoColors(), raster =F)


####
#### Fig. S3C  - Violin plots of quality control metrics for WIBR3 DAN dataset.
####
dittoPlot(all, "nFeature_RNA", group.by = "subclone_ID", vlnplot.width = 0.5, plots = c("vlnplot","boxplot"), boxplot.show.outliers = T, vlnplot.lineweight = 0.25, boxplot.lineweight = 0.25)
dittoPlot(all, "nCount_RNA", group.by = "subclone_ID", vlnplot.width = 0.5, plots = c("vlnplot","boxplot"), boxplot.show.outliers = T, vlnplot.lineweight = 0.25, boxplot.lineweight = 0.25)


####
#### Fig. S3D  - FeaturePlots of key marker genes in the WIBR3 DAN dataset.
####
FeaturePlot(all, features = "TH", reduction = "umap", order = TRUE, pt.size = 2)
FeaturePlot(all, features = "KCNJ6", reduction = "umap", order = TRUE, pt.size = 2)
FeaturePlot(all, features = "NR4A2", reduction = "umap", order = TRUE, pt.size = 2)
FeaturePlot(all, features = "CALB1", reduction = "umap", order = TRUE, pt.size = 2)
FeaturePlot(all, features = "SOX6", reduction = "umap", order = TRUE, pt.size = 2)
FeaturePlot(all, features = "FOXA2", reduction = "umap", order = TRUE, pt.size = 2)


###
### Fig. S4A-i - UMAP of WIBR3 DAN cells, with FOUNDIN-PD cell type labels (applied by SingleR with FOUNDIN-PD cell types used as a reference).
###
dittoDimPlot(all, "FOUNDIN_PD_assignment_singleR",
             colors = c(10,1,6,4,3,11,5,8,2,9,7), 
             do.label = T, 
             labels.repel = T, 
             labels.size = 5, 
             size = 1.5)

###
### Fig. S4A-ii - Stacked barplots showing the breakdown of cell types in both the WIBR3 DAN and FOUNDIN-PD datasets.
###
### First for WIBR3 DAN dataset.
Idents(all) <- "FOUNDIN_PD_assignment_singleR"
pt <- table(Idents(all), all$FOUNDIN_PD_assignment_singleR)
pt <- as.data.frame(pt)
pt <- filter(pt, Freq > 0)
desired_order <- c("iDA1", "iDA2", "iDA3", "iDA4", "lProg1", "lProg2", "eProg1", "eProg2", "PFPP", "Ependymal", "NE")
pt$Var1 <- factor(pt$Var1, levels = desired_order)
pt <- arrange(pt, factor(Var1, levels = c("iDA1","iDA2","iDA3","iDA4","lProg1","lProg2","eProg1","eProg2","PFPP","Ependymal","NE")))
ggplot(pt, aes(x = 1, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) + geom_col(position = "fill", width = 0.5) +
    xlab("Sample") + ylab("Proportion") +
    scale_fill_manual(values=colors) + theme(legend.title = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
### Repeat for FOUNDIN-PD dataset.
Idents(foundin) <- "CellType"
pt <- table(Idents(foundin), foundin$CellType)
pt <- as.data.frame(pt)
pt <- filter(pt, Freq > 0)
desired_order <- c("iDA1", "iDA2", "iDA3", "iDA4", "lProg1", "lProg2", "eProg1", "eProg2", "PFPP", "Ependymal", "NE")
pt$Var1 <- factor(pt$Var1, levels = desired_order)
pt <- arrange(pt, factor(Var1, levels = c("iDA1","iDA2","iDA3","iDA4","lProg1","lProg2","eProg1","eProg2","PFPP","Ependymal","NE")))
ggplot(pt, aes(x = 1, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) + geom_col(position = "fill", width = 0.5) +
    xlab("Sample") + ylab("Proportion") +
    scale_fill_manual(values=colors) + theme(legend.title = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


###
### Fig. S4A-iii - Default umap of the FOUNDIN-PD dataset, with cell type labels applied.
###
dittoDimPlot(foundin, "CellType", 
    reduction.use = "umap", 
    do.label = T, 
    labels.repel =T,  
    do.raster = TRUE)


###
### Fig. S4B - Heatmap of SingleR cell assignment scores, comparing WIBR3 DAN cells to the full FOUNDIN-PD reference dataset.
###
FOUNDIN_PD_assignment_singleR <- readRDS("SingleR_FOUNDIN-PD_assignments.rds")
plotScoreHeatmap(FOUNDIN_PD_assignment_singleR, 
                 show.labels = T, 
                 cluster_cols = F, 
                 show.pruned = F)


###
### Fig. S4C - WIBR3 DAN cells, labeled with FOUNDIN-PD cell type labels, presented in the UMAP space from the FOUNDIN-PD dataset.
###
dittoDimPlot(all, "FOUNDIN_PD_assignment_singleR", 
    reduction.use = "ref.umap", 
    colors = c(10,1,6,4,3,11,5,8,2,9,7), 
    do.label = T, 
    labels.repel = T, 
    labels.size = 5, 
    size = 1.5)


###
### Fig. S4D - Heatmap of top 10 marker genes that define the FOUNDIN-PD cell type labels that were applied to the WIBR3 DAN dataset by SingleR.
###
markers(all, "all_WIBR3_DAN_SingleR_Assignments_MAST", "FOUNDIN_PD_assignment_singleR", 0.3, 0.5, "MAST")
m <- read.csv("all_WIBR3_DAN_SingleR_Assignments_MAST_markers_by_FOUNDIN_PD_assignment_singleR.csv")
m %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(all, features = top10$gene, group.colors = c("#F0E442","#009E73","#007756","#0072B2",
                                                           "#666666","#56B4E9","#E69F00","#D55E00",
                                                           "#CC79A7","#1C91D4","#AD7700"),raster =F)