#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name="HyprColoc"

# Here load needed system tools (Java 1.8 is required, one of singularity or anaconda - python 2.7 are needed,
# depending on the method for dependancy management)

module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4

nextflow_path=[Full path to the folder where Nextflow is installed]

${nextflow_path}/nextflow run RunHyprColocOnGWAS.nf \
--gwas '[Path to GWAS summary statistics file]' \
--window '[Window around lead SNP to search colocalization]' \
--Pthresh '[P-value threshold for lead SNP]' \
--OutputDir '[Output folder]' \
-with-report HyprColoc.html \
-resume \
-profile singularity,slurm
