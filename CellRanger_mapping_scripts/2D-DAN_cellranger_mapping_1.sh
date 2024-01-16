#!/bin/bash
#SBATCH --time 0-24 # days-minutes
#SBATCH --job-name=2D-DAN_CRmulti-1 # Job name
#SBATCH --account=fc_hockemeyer
#SBATCH --nodes=1
#SBATCH --ntasks=20 # Number of cores
#SBATCH --mem=100000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=savio2_bigmem # Partition to submit to
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jessedunnack@berkeley.edu
#SBATCH --output=2D-DAN_CRmulti-1_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=2D-DAN_CRmulti-1_%A_%a.err # File to which STDERR will be written


##making sure the path is set
export PATH=/global/scratch/jessedunnack/cellranger-7.0.1:$PATH

cellranger multi --id=Part_1 --csv=/global/scratch/users/jessedunnack/2D_DAN_scRNASeq/2D-DAN_sample_info_1.csv
cellranger multi --id=Part_2 --csv=/global/scratch/users/jessedunnack/2D_DAN_scRNASeq/2D-DAN_sample_info_2.csv
cellranger multi --id=Part_3 --csv=/global/scratch/users/jessedunnack/2D_DAN_scRNASeq/2D-DAN_sample_info_3.csv
