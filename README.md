# iSCORE-PD_DAN_scRNASeq

**iSCORE-PD: an isogenic stem cell collection to research Parkinson’s Disease**

**Oriol Busquets**, **Hanqin Li**, **Khaja Mohieddin Syed**, **Pilar Alvarez Jerez**, **Jesse Dunnack**, Riana Lo Bu, Yogendra Verma, Gabriella R. Pangilinan, Annika Martin, Jannes Straub, YuXin Du, Vivien M. Simon, Steven Poser, Zipporiah Bush, Jessica Diaz, Atehsa Sahagun, Jianpu Gao, Dena G. Hernandez, Kristin S. Levine, Ezgi O. Booth, Helen S. Bateup, Donald C. Rio, Dirk Hockemeyer, Cornelis Blauwendraat, Frank Soldner
bioRxiv 2024.02.12.579917; doi: https://doi.org/10.1101/2024.02.12.579917

**WORKFLOW OVERVIEW: **

In the iSCORE-PD manuscript above, this collection of scripts is used to analyze scRNASeq data from hESC-derived dopaminergic neurons.

    A) CellRanger:
      1) Map raw reads + generate count matrices with 10x Genomics' CellRanger suite.
    
    B) Seurat in R:
      1) Initialize Seurat objects for each sample, perform QC filtering on each Seurat object.
      2) Integrate the 3 subclones to generate an aggregate dataset of human midbrain DANs.
      3) Evaluate cells against FOUNDIN-PD reference dataset.
      4) Visualize results, generate necessary plots.

[![DOI](https://zenodo.org/badge/741244123.svg)](https://zenodo.org/doi/10.5281/zenodo.10718769)
