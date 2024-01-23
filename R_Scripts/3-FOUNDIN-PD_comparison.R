## Run on HPC cluster to speed up analysis
## Used conda environment "seuratv4"
## See "seuratv4_environment.yml" for more info.

library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(scuttle)
library(scran)
library(dplyr)

options(future.globals.maxSize = 10000 * 1024^2)

reclust <- function(x, vars_to_regress){
  x <- SCTransform(x, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = vars_to_regress)
  x <- RunPCA(x, npcs = 100)
  x <- FindNeighbors(x, dims = 1:100)
  x <- FindClusters(x)
  x <- RunUMAP(x, dims = 1:100, return.model = T)
  return(x)
}


###
### FOUNDIN-PD comparison using Seurat
###

all <- readRDS("all_WIBR3_DAN.rds") 
foundin <- readRDS("../FOUNDIN-PD/iPSCsDopaALL_integratedAfterBroadCellType.RDS")
foundin <- RunUMAP(foundin, dims = 1:30, return.model=TRUE)
all.anchors <- FindTransferAnchors(reference = foundin, query = all, dims = 1:30,
                                    reference.reduction = "pca", normalization.method = "LogNormalize", 
                                    query.assay = "RNA", reference.assay = "RNA", features = intersect(rownames(all),rownames(foundin)))
predictions <- TransferData(anchorset = all.anchors, refdata = foundin$CellType, dims = 1:30)
saveRDS(predictions, "Seurat_FOUNDIN-PD_assignments.rds")

all <- AddMetaData(all, metadata = predictions)
all <- MapQuery(anchorset = all.anchors, reference = foundin, query = all,
                refdata = list(celltype = "CellType"), reference.reduction = "pca", reduction.model = "umap")

# Append FOUNDIN-PD cell assignments done by Seurat to our dataset.
predictions <- readRDS("Seurat_FOUNDIN-PD_assignments.rds")
all$FOUNDIN_PD_assignment_seurat <- predictions$predicted.id
all$FOUNDIN_PD_assignment_score_seurat <- predictions$prediction.score.max
saveRDS(all, "all_WIBR3_DAN.rds")

###
### FOUNDIN-PD comparison using SingleR
###

all <- readRDS("all_WIBR3_DAN.rds")
all.sce <- logNormCounts(as.SingleCellExperiment(all))
foundin <- readRDS("../FOUNDIN-PD/iPSCsDopaALL_integratedAfterBroadCellType.RDS")
foundin <- logNormCounts(as.SingleCellExperiment(foundin))

pred_wilcox <- SingleR(test= all.sce,
                        ref = foundin,
                        labels = foundin$CellType,
                        de.n = 100,
                        de.method = "wilcox")
saveRDS(pred_wilcox, "SingleR_FOUNDIN-PD_assignments.rds")


# Append FOUNDIN-PD cell assignments done by SingleR to our dataset.
pred_wilcox <- readRDS("SingleR_FOUNDIN-PD_assignments.rds") 
all$FOUNDIN_PD_assignment_singleR <- as.data.frame(pred_wilcox)$pruned.labels

meta_data_table <- all@meta.data
singleR_table <- as.data.frame(FOUNDIN_PD_assignment_singleR)
for (row_index in 1:nrow(meta_data_table)) {
    x <- meta_data_table$FOUNDIN_PD_assignment_singleR[row_index]
    if (!is.na(x)) {
        score_column_name <- paste0("scores.", x)
        score_entry <- singleR_table[row_index, score_column_name]
        if (length(score_entry) == 0) {
            score_entry <- NA 
        }
                meta_data_table$FOUNDIN_PD_assignment_score_singleR[row_index] <- score_entry
    } else {
        meta_data_table$FOUNDIN_PD_assignment_score_singleR[row_index] <- NA
    }
}
all@meta.data <- meta_data_table

## Subset to remove any cells that are "pruned" by SingleR during the FOUNDIN-PD comparison
rows_with_na <- which(is.na(all$FOUNDIN_PD_assignment_singleR))
all <- all[, -rows_with_na]

all <- reclust(all, c("percent.mt","nCount_CMO"))
all$seurat_clusters_fine <- all$seurat_clusters

all <- FindClusters(all, resolution = 0.18)
all$seurat_clusters_coarse <- all$seurat_clusters

#Create another column that normalizes all SingleR assignment scores to 1.
x <- all$FOUNDIN_PD_assignment_score_singleR
normalized = (x-min(x))/(max(x)-min(x))
all@meta.data$FOUNDIN_PD_assignment_score_singleR_normalized <- normalized


# Clean up metadata of final Seurat object.
all$seurat_clusters <- NULL
all$SCT_snn_res.0.18 <- NULL
all$SCT_snn_res.0.8 <- NULL
all@meta.data <- all@meta.data %>%
  select(c("experiment_ID","genetic_background","genotype_ID","subclone_ID","age",
  "assay_date","nCount_RNA","nFeature_RNA","nCount_CMO","nFeature_CMO","nCount_SCT",
  "nFeature_SCT","percent.mt","percent.ribo","S.Score","G2M.Score","Phase",
  "CR_Multiplet","CR_Assignment","CR_Assignment_Probability","seurat_clusters_fine",
  "seurat_clusters_coarse","FOUNDIN_PD_assignment_seurat","FOUNDIN_PD_assignment_score_seurat",
  "FOUNDIN_PD_assignment_singleR","FOUNDIN_PD_assignment_score_singleR",
  "FOUNDIN_PD_assignment_score_singleR_normalized"))

saveRDS(all, "all_WIBR3_DAN.rds")