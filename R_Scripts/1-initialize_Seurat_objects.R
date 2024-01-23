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

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
setwd("E:/ASAP/scRNASeq/2D_DAN/all_WT_analysis/iSCORE-PD/")

add_cellranger_demux_results <- function(seurat_obj, name){
    assignment_table <- read.csv(paste0("./cellRanger_outputs/",name,"/",name,"_assignment_confidence_table.csv"))
    rownames(assignment_table) <- assignment_table$Barcode
    assignment_table <- subset(assignment_table, select = c("Multiplet", "Assignment", "Assignment_Probability"))
    colnames(assignment_table) <- paste0('CR_', colnames(assignment_table))
    seurat_obj <- AddMetaData(seurat_obj, metadata = assignment_table)
    return(seurat_obj)
}
add_metadata <- function(seurat_obj, name, meta_csv_path) {
    meta <- read.csv(meta_csv_path, header = TRUE)
    matching_rows <- grepl(name, meta$experiment_ID)
    submeta <- meta[matching_rows, ]
    experiment_cells <- grepl(name, seurat_obj$experiment_ID)
    matching_CMO_ID <- seurat_obj$CR_Assignment[experiment_cells]
    matching_rows <- match(matching_CMO_ID, submeta$CMO_ID)
    seurat_obj@meta.data[experiment_cells, "subclone_ID"] <- submeta$subclone_ID[matching_rows]
    seurat_obj@meta.data[experiment_cells, "genotype_ID"] <- submeta$genotype_ID[matching_rows]
    seurat_obj@meta.data[experiment_cells, "genetic_background"] <- submeta$genetic_background[matching_rows]
    seurat_obj@meta.data[experiment_cells, "age"] <- submeta$age[matching_rows]
    seurat_obj@meta.data[experiment_cells, "assay_date"] <- submeta$assay_date[matching_rows]
    return(seurat_obj)
}
initialize_seurat_object <- function(name) {
    raw_counts <- Read10X(paste0("./count_matrices/",name,"/"))
    assignment_table <- read.csv(paste0("./cellRanger_outputs/",name,"/",name,"_assignment_confidence_table.csv"))
    barcodes_to_keep <- assignment_table$Barcode
    x <- CreateSeuratObject(counts = raw_counts$`Gene Expression`)
    x[['CMO']] = CreateAssayObject(counts = raw_counts$`Multiplexing Capture`)
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
    x[["percent.ribo"]] <- PercentageFeatureSet(x, pattern = "^RP[SL]")
    x[["experiment_ID"]] <- name
    x <- subset(x, cells = barcodes_to_keep)
    x <- CellCycleScoring(object = x, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
    x <- add_cellranger_demux_results(x, name)
    x <- add_metadata(x, name, "./iSCORE-PD_scRNASeq_METADATA.csv")
    x <- subset(x, CR_Assignment_Probability >= 0.9)
    x <- subset(x, genotype_ID == "WT")
    x <- subset(x, nFeature_RNA >= 1500 & nCount_RNA < 30000 & percent.mt <= 10)
    x <- SCTransform(x, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = "percent.mt")
    x <- RunPCA(x, npcs = 100)
    x <- FindNeighbors(x, dims = 1:100)
    x <- FindClusters(x)
    x <- RunUMAP(x, dims = 1:100)
    return(x)
}

WIBR3_S1_DAN <- initialize_seurat_object("DAN_01")
WIBR3_S2_DAN <- initialize_seurat_object("DAN_02")
WIBR3_S3_DAN <- initialize_seurat_object("DAN_03")

# Rearrange metadata columns to tidy things up.
WIBR3_S1_DAN@meta.data <- WIBR3_S1_DAN@meta.data %>%
  select(c("experiment_ID","genetic_background","genotype_ID","subclone_ID","age",
  "assay_date","nCount_RNA","nFeature_RNA","nCount_CMO","nFeature_CMO","nCount_SCT",
  "nFeature_SCT","percent.mt","percent.ribo","S.Score","G2M.Score","Phase",
  "CR_Multiplet","CR_Assignment","CR_Assignment_Probability","seurat_clusters"))

WIBR3_S2_DAN@meta.data <- WIBR3_S2_DAN@meta.data %>%
  select(c("experiment_ID","genetic_background","genotype_ID","subclone_ID","age",
  "assay_date","nCount_RNA","nFeature_RNA","nCount_CMO","nFeature_CMO","nCount_SCT",
  "nFeature_SCT","percent.mt","percent.ribo","S.Score","G2M.Score","Phase",
  "CR_Multiplet","CR_Assignment","CR_Assignment_Probability","seurat_clusters"))

WIBR3_S3_DAN@meta.data <- WIBR3_S3_DAN@meta.data %>%
  select(c("experiment_ID","genetic_background","genotype_ID","subclone_ID","age",
  "assay_date","nCount_RNA","nFeature_RNA","nCount_CMO","nFeature_CMO","nCount_SCT",
  "nFeature_SCT","percent.mt","percent.ribo","S.Score","G2M.Score","Phase",
  "CR_Multiplet","CR_Assignment","CR_Assignment_Probability","seurat_clusters"))


saveRDS(WIBR3_S1_DAN, "WIBR3_S1_DAN.rds")
saveRDS(WIBR3_S2_DAN, "WIBR3_S2_DAN.rds")
saveRDS(WIBR3_S3_DAN, "WIBR3_S3_DAN.rds")