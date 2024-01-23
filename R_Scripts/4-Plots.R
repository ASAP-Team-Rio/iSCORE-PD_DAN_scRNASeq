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
markers <- function(x, phrase, var, min_pct, lfc_cutoff) {
     if (is.null(var)) {
         markers <- FindAllMarkers(x, only.pos = FALSE, min.pct = min_pct, logfc.threshold = lfc_cutoff) 
         markers %>%
             group_by(cluster) %>%
             top_n(n = 10, wt = avg_log2FC) -> top10
         p <- DoHeatmap(x, features = top10$gene) + NoLegend()
         ggsave(filename = paste0(phrase, "_MarkerHeatMap.png"), plot = p, device = "png", width = 30, height = 25)
         write.csv(markers, file = paste0(phrase, "_markers.csv"))
         write.csv(top10, file = paste0(phrase, "_top10markers.csv"))
     } else {
         x <- SetIdent(x, value = var)
         markers <- FindAllMarkers(x, only.pos = FALSE, min.pct = min_pct, logfc.threshold = lfc_cutoff) 
         markers %>%
             group_by(cluster) %>%
             top_n(n = 10, wt = avg_log2FC) -> top10
         p <- DoHeatmap(x, features = top10$gene) + NoLegend()
         ggsave(filename = paste0(phrase, "_MarkerHeatMap.png"), plot = p, device = "png", width = 30, height = 25)
         write.csv(markers, file = paste0(phrase, "_markers_by_", var, ".csv"))
         write.csv(top10, file = paste0(phrase, "_top10markers_by_", var, ".csv")) 
     }
}

### PSEUDOBULK SPEARMAN CODE FIG 2E
mtx <- pseudoBulk_spearman(all, foundin, "seurat_clusters_coarse", "CellType", "RNA", n_var_features = 3000)
saveRDS(mtx, "all_v_FOUNDIN_FULL_Spearman_matrix_RNA3000.rds")
mtx <- readRDS("all_v_FOUNDIN_Spearman_matrix_RNA3000.rds")
pheatmap(mtx2,
    main = "Heatmap of Spearman correlations \n WIBR3 vs. FOUNDIN-PD DA Neurons\n (RNA Assay, pseudobulked clusters) \n Input = 3000 variable genes", 
    fontsize = 12, 
    clustering_method = "ward.D2", 
    cutree_rows = 4, 
    cutree_cols = 3, 
    treeheight_row = 30, 
    treeheight_col = 30)


#Identifying marker genes that define identity classes
markers(all, "all_WIBR3_DAN_coarse_clusters", "seurat_clusters_coarse", 0.3, 0.5)
markers(all, "all_WIBR3_DAN_fine_clusters", "seurat_clusters_fine", 0.3, 0.5)
markers(all, "all_WIBR3_DAN_SingleR_Assignments", "FOUNDIN_PD_assignment_singleR", 0.3, 0.5)
markers(all, "all_WIBR3_DAN_Seurat_Assignments", "FOUNDIN_PD_assignment_seurat", 0.3, 0.5)

#Plotting the heatmap
m <- read.csv("all_WIBR3_DAN_markers_by_seurat_clusters.csv")
m %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
p <- DoHeatmap(all, features = top10$gene) + NoLegend()
p
ggsave(filename = "all_WIBR3_DAN_Coarse_Clusters_Top10_MarkerHeatMap.png", plot = p, device = "png", width = 30, height = 25)



###
#DOT PLOT CODE
###
goi <- c("GAP43", "MAPT","STMN2","TUBB3","TAC1","NEFL","SNCA","NCAM2",
         "ZCCHC12","ATP1A3","SNAP25","VGF","NKAIN2","CNTNAP2","CTNNA2",
         "SYT1","NRXN1","NRG3","SGCZ","NRG1","DCC","CNTNAP5","KCNJ6","NRXN3","TH",
         "SLC18A2","SLC18A1","TPH1","CALCA","CGA","CHGA","DLK1",
         "VCAN","LGALS1","HMGB2","CENPF","UBE2C","TOP2A","MKI67",
         "HES1","CALB1","SLIT2","VIM","SPARC","FOXP2","CORIN","SOX6")
goi_2 <- c("MAPT","STMN2","GAP43","SNAP25","SNCA","CNTNAP2","LMX1A","LMX1B","KCNJ6",
           "TH", "SLC18A1","SLC18A2","TPH1","CALCA","CGA","CHGA","DLK1",
           "VCAN","LGALS1","UBE2C","TOP2A","MKI67","HES1","CALB1","SLIT2",
           "VIM","SPARC","FOXP2","CORIN","SOX6","FOXA2","NR4A2")
dittoDotPlot(all, goi_2, group.by = "seurat_clusters", y.reorder =c (1,2,5,4,6,3), size = 8)

#SEURAT DOTPLOT VERSION
all@active.ident <- factor(all@active.ident, levels=c("0","1","4","6","3","5","2"))
goi_2 <- c("STMN2","GAP43","MAPT","SNAP25","SNCA","CNTNAP2","LMX1A","LMX1B","KCNJ6",
           "TH","NR4A2","FOXA2","SLC18A1","SLC18A2","TPH1","CALCA","CGA","CHGA","DLK1",
           "VCAN","LGALS1","UBE2C","TOP2A","MKI67","HES1","CALB1","SLIT2",
           "VIM","SPARC","FOXP2","CORIN","SOX6")
DotPlot(all, features = goi_2, assay = 'SCT', dot.scale = 8, scale.by = "size")+scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "Oranges")) +RotatedAxis()


all@active.ident <- factor(all@active.ident, levels=c("0","1","4","6","3","5","2"))
goi_2 <- c("TH","NR4A2","KCNJ6","LMX1B","LMX1A","STMN2",
           "MAP2","GAP43","SNAP25","SNCA","SYN1","ROBO1",
           "FOXA2","SLC18A1","SLC18A2","TPH1","MKI67",
           "HES1","SLIT2","FOXP2","CORIN","SOX6","CALB1")
DotPlot(all, features = goi_2, assay = 'SCT', dot.scale = 8, scale.by = "size")+scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "Oranges")) +RotatedAxis()



#FINAL UMAP PLOT CODE:

a <- dittoDimPlot(all, "seurat_clusters", reduction.use = "ref.umap", main = "WIBR3 DAN", sub = "By Seurat Clustering", do.contour = T, contour.linetype = 6)
b <- dittoDimPlot(all, all$SingleR_wilcox100_pruned_labels, reduction.use = "ref.umap", main = "WIBR3 DAN", sub = "By SingleR cell type assignments", do.contour = T, contour.linetype = 6)
c <- dittoDimPlot(all, all$Seurat_label_transfer_celltype, reduction.use = "ref.umap", main = "WIBR3 DAN", sub = "By Seurat Label Transfer", do.contour = T, contour.linetype = 6)
wrap_plots(a,a,b,c, ncol = 2, nrow = 2)



x <- dittoDimPlot(all, "seurat_clusters", reduction.use = "umap", main = "WIBR3 DAN", sub = "By Seurat Clustering", do.contour = T, contour.linetype = 3)
y <- dittoDimPlot(all, all$SingleR_wilcox100_pruned_labels, reduction.use = "umap", main = "WIBR3 DAN", sub = "By SingleR cell type assignments", do.contour = T, contour.linetype = 3)
z <- dittoDimPlot(all, all$Seurat_label_transfer_celltype, reduction.use = "umap", main = "WIBR3 DAN", sub = "By Seurat Label Transfer", do.contour = T, contour.linetype = 3)
wrap_plots(x,x,y,z, ncol = 2, nrow = 2)



dittoDimPlot(all, "seurat_clusters_coarse", do.label = T, labels.size = 12, labels.repel = T)
#### UMAP SPLIT BY SUBCLONE
dittoDimPlot(all, "seurat_clusters_coarse", split.by = "subclone_ID", do.label = T, labels.size = 12, labels.repel = T)

##Violin Plot code ffor supplement:
dittoPlot(x, "nFeature_RNA", group.by = "subclone_ID", vlnplot.width = 0.5, plots = c("vlnplot","boxplot"), boxplot.show.outliers = T, vlnplot.lineweight = 0.25, boxplot.lineweight = 0.25)
dittoPlot(x, "nCount_RNA", group.by = "subclone_ID", vlnplot.width = 0.5, plots = c("vlnplot","boxplot"), boxplot.show.outliers = T, vlnplot.lineweight = 0.25, boxplot.lineweight = 0.25)


### DimPlots that show FOUNDIN-PD cell assignment scores
dittoDimPlot(all, "FOUNDIN_PD_assignment_score_seurat", size = 1.5)
dittoDimPlot(all, "FOUNDIN_PD_assignment_score_singleR", size = 1.5)
dittoDimPlot(all, "FOUNDIN_PD_assignment_score_singleR_normalized", size = 1.5)




plots(all, "all_WIBR3_DAN")
markers(all, "seurat_clusters_coarse", "all_WIBR3_DAN_coarse_clusters", 0.3, 0.5)
markers(all, "seurat_clusters_fine", "all_WIBR3_DAN_fine_clusters", 0.3, 0.5)
markers(all, "FOUNDIN_PD_assignment_seurat", "all_WIBR3_DAN_FOUNDIN-PD_labels_seurat", 0.3, 0.5)
markers(all, "FOUNDIN_PD_assignment_singleR", "all_WIBR3_DAN_FOUNDIN-PD_labels_singleR", 0.3, 0.5)




#Figure Code - SingleR score heatmap
FOUNDIN_PD_assignment_singleR <- readRDS("SingleR_FOUNDIN-PD_assignments.rds")

plotScoreHeatmap(FOUNDIN_PD_assignment_singleR, 
                 show.labels = T, 
                 cluster_cols = F, 
                 show.pruned = F)




c("Ependymal","eProg1","eProg2","iDA1","iDA2","iDA3","iDA4","lProg1","lProg2","NE","PFPP")
