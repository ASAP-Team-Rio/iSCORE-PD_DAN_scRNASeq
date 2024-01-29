# iSCORE-PD_DAN_scRNASeq

iSCORE-PD: An isogenic stem cell collection to research Parkinson’s Disease

The following is a codebase covering the analysis of scRNASeq data of human midbrain dopaminergic neurons differentiated in vitro from hESCs.

OVERVIEW:
A) CellRanger:
  1) Map raw reads + generate count matrices with 10x Genomics' CellRanger suite.

B) Seurat in R:
  1) Initialize Seurat objects for each sample, perform QC filtering on each Seurat object.
  2) Integrate the 3 subclones to generate an aggregate dataset of human midbrain DANs.
  3) Evaluate cells against FOUNDIN-PD reference dataset.
  4) Visualize results, generate necessary plots.
