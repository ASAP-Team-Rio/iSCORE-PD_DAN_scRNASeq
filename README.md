# iSCORE-PD_DAN_scRNASeq

iSCORE-PD: An isogenic stem cell collection to research Parkinson’s Disease

The following is a codebase covering the analaysis of scRNASeq data of human midbrain dopaminergic neurons differentiated in vitro from hESCs.

OVERVIEW:
1) Map raw reads + generate count matrices with 10x Genomics' CellRanger suite.
2) Initialize Seurat objects for each sample.
3) Perform cell demux analysis with cellHashR.
4) Perform QC filtering on individual Seurat objects.
5) Integrate the 3 subclones to generate an aggregate dataset of human midbrain DANs.
6) Evaluate cells against FOUNDIN-PD reference dataset.
