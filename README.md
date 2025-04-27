# iSCORE-PD_DAN_scRNASeq

**iSCORE-PD: an isogenic stem cell collection to research Parkinson’s Disease**

**Oriol Busquets**, **Hanqin Li**, **Khaja Mohieddin Syed**, **Pilar Alvarez Jerez**, **Jesse Dunnack**, Riana Lo Bu, Yogendra Verma, Gabriella R. Pangilinan, Annika Martin, Jannes Straub, YuXin Du, Vivien M. Simon, Steven Poser, Zipporiah Bush, Jessica Diaz, Atehsa Sahagun, Jianpu Gao, Dena G. Hernandez, Kristin S. Levine, Ezgi O. Booth, Helen S. Bateup, Donald C. Rio, Dirk Hockemeyer, Cornelis Blauwendraat, Frank Soldner
bioRxiv 2024.02.12.579917; doi: https://doi.org/10.1101/2024.02.12.579917
[![DOI](https://zenodo.org/badge/741244123.svg)](https://zenodo.org/doi/10.5281/zenodo.10718769)




**WORKFLOW OVERVIEW:**

In the iSCORE-PD manuscript above, this collection of scripts is used to analyze scRNASeq data from hESC-derived dopaminergic neurons. CellRanger was run on Linux in an HPC environment, then count matrices were processed and analyzed in R on Windows 11.

A) CellRanger:
   1) Map raw reads, demultiplex samples and generate count matrices with 10x Genomics' CellRanger suite.
      NOTE - Cell multiplexing oligos (CMO) associated with samples are specified in the "SampleInfo" CSVs.
    
B) Seurat in R:
   1) Initialize Seurat objects for each sample, perform QC filtering on each Seurat object.
   2) Integrate the 3 subclones to generate an aggregate dataset of human midbrain DANs.
   3) Evaluate cells against FOUNDIN-PD reference dataset.
   4) Visualize results, generate necessary plots.


#### Quick start
We provide a Conda environment file to set up all required dependencies. In Linux:

```bash
git clone https://github.com/ASAP-Team-Rio/iSCORE-PD_DAN_scRNASeq.git
cd iSCORE-PD_DAN_scRNASeq/

conda env create -f R_Scripts/seuratv4_environment.yml
conda activate iSCORE-PD-scRNAseq
```

### Data Download
The RDS files containing processed Seurat objects can be accessed at [URL_TO_DATA].

## Instructions for Use

### Full Analysis Workflow
1. Process raw data:
   ```R
   Rscript R_Scripts/01_preprocessing.R
   ```

2. Perform integration and clustering:
   ```R
   Rscript R_Scripts/02_integration_clustering.R
   ```

3. Generate figures:
   ```R
   Rscript R_Scripts/03_visualization.R
   ```

## System Requirements

### Hardware Requirements
- RAM: 16GB minimum, 32GB+ recommended for larger datasets
- CPU: Multi-core processor (8+ cores recommended)

### Operating Systems
- Tested on:
  - Linux (Rocky Linux release 8.10 - Green Obsidian)
  - Windows 11 - x64 (build 22631)

### Dependencies
- This code requires R (version 4.1.0 or higher) with various packages, specified in "/R_Scripts/R_sessioninfo.txt"
- Cell Ranger (10x Genomics) - v8.0.1




## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
For questions about the code, please contact Jesse Dunnack - jessedunnack@gmail.com / jessedunnack@berkeley.edu
