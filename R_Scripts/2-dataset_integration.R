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

WIBR3_S1_DAN <- readRDS("WIBR3_S1_DAN.rds")
WIBR3_S2_DAN <- readRDS("WIBR3_S2_DAN.rds")
WIBR3_S3_DAN <- readRDS("WIBR3_S3_DAN.rds")

reclust <- function(x, vars_to_regress){
  x <- SCTransform(x, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = vars_to_regress)
  x <- RunPCA(x, npcs = 100)
  x <- FindNeighbors(x, dims = 1:100)
  x <- FindClusters(x)
  x <- RunUMAP(x, dims = 1:100, return.model = T)
  return(x)
}

all <- list(WIBR3_S1_DAN, WIBR3_S2_DAN, WIBR3_S3_DAN)
all_features <- SelectIntegrationFeatures(object.list = all, nfeatures = 3000)
all <- PrepSCTIntegration(object.list = all, anchor.features = all_features)
all_anchors <- FindIntegrationAnchors(object.list = all, normalization.method = "SCT", anchor.features = all_features)
all <- IntegrateData(anchorset = all_anchors, normalization.method = "SCT")
all <- reclust(all, c("percent.mt","nCount_CMO"))

saveRDS(all, "all_WIBR3_DAN.rds")