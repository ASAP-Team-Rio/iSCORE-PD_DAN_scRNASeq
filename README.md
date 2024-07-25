# iSCORE-PD_DAN_scRNASeq

iSCORE-PD: An isogenic stem cell collection to research Parkinson’s Disease

In the iSCORE-PD manuscript above, this collection of scripts is used to analyze scRNASeq data from hESC-derived dopaminergic neurons.

WORKFLOW OVERVIEW: 

A) CellRanger:
  1) Map raw reads + generate count matrices with 10x Genomics' CellRanger suite.

B) Seurat in R:
  1) Initialize Seurat objects for each sample, perform QC filtering on each Seurat object.
  2) Integrate the 3 subclones to generate an aggregate dataset of human midbrain DANs.
  3) Evaluate cells against FOUNDIN-PD reference dataset.
  4) Visualize results, generate necessary plots.

[![DOI](https://zenodo.org/badge/741244123.svg)](https://zenodo.org/doi/10.5281/zenodo.10718769)
