# Pipeline for running colocalisation

This pipeline takes GWAS summary statistics and runs HyprColoc analyses for the datasets in eQTL Catalogue.

## Word of caution

By default, this pipeline extracts all associated loci from your input GWAS and tests each of those loci against 433 eQTL Catalogue datasets. Each pipeline step is a separate job which is submitted into the scheduling system. Therefore, if you want to analyse GWAS with 100 associated loci, you will submit 1 + 100 + 433 * 100 * 2 = **86,701** jobs into the HPC queue (max 400 at the same time). Although each individual job requires small resources and completes relatively quickly, it is still considerable amount of calculations. Therefore it is recommended:

1. Use default settings only if your GWAS has not too many associated loci.
2. If you have very polygenic GWAS (e.g. hundreds of associated loci), consider limiting your analysis to the smaller subset of most interesting loci you want to analyse in-depth.
3. Or limit the eQTL Catalogue datasets you include to the biologically relevant tissues/cell types.
4. Or include only certain most relevant transcriptional features from eQTL Catalogue (e.g. gene expression QTLs but not splicing, exon, txrev. This would limit the number of jobs to 25% from the full analysis).
5. Use some hypothesis-based approach to bring down the number of loci and tissues you want to analyse.

## Instructions

### Requirements

- You need to have Java 8 installed in your HPC environment.
- You need to have Singularity installed in your HPC environment.
- You need SLURM as scheduler.
- You need to download and install Nextflow executable (https://www.nextflow.io/docs/latest/getstarted.html#installation).

### Setup of the analysis

In order to keep analysis environment clean and tidy, I recommend to run all Nextflow pipelines in the following structure.

1. Make analysis folder. E.g. `HeightGwasHyprcoloc`.
2. Inside this folder make needed folders for inputs and outputs. E.g. `input` and `output`. In case help files are needed, make separate `helpfiles`. You can make `tools` folder and move your Nextflow executable into this folder. Additional folders can be made, according to need. Move all the inputs and help files into those folders.
3. Download pipeline from gitlab. Best way: use `git clone`. As a result you should have separate folder for the nextflow pipeline.
4. In the downloaded folder is script template. Adjust the inputs/outputs in this script template according to your inputs. I also rename the script after modification.
5. Submit the modified script from pipeline folder.

### Input

By default expects GWAS summary statistics to be in GWAMA output format and expects that these are in hg38 (format of eQTL catalogue). However it is possible to specify the input column names in `conf/params.config` or with argument flags.

### Settings

**Mandatory arguments:**

`--gwas`          GWAS summary statistics file.

**Optional arguments:**

`--GwasType`                Whether your GWAS trait is binary (case-control) or continuous. One of two: `cc` or `cont`. Defaults to `cc`

`--OutputDir`               Directory for colocalisation results and SNP PIPs from each GWAS-eQTL colocalisation (colocalisation posterior > 0.5).

`--window`                  Window around each lead SNP to test the colocalization. Defaults to 1000000 (+/-1Mb from lead SNP). It is not recommended to use larger window.

`--Pthresh`                 P-value threshold to identify lead SNPs. Defaults to 5e-8.

`--MafThresh`               MAF threshold to filter the input GWAS data. Defaults to 0.001.

`--PosteriorThreshold`      Posterior probability to declare colocalization. Defaults to 0.8.

`--CsThreshold`             Threshold for credible set calculation. Defaults to 0.95 (95% credible sets).

`--OutputCsPip`             Whether to write out PIPs for each colocalising locus: TRUE/FALSE, defaults to FALSE.

`--eQtlFile`                Local file with eQTL Catalogue links. If this is specified, only those files are used to run the analysis. It is useful for limiting the analysis to the datasets of interest and limiting the number of jobs submitted to the queue. Secondly, this file can be used to run analysis on the local data folder (not to download from ftp). In the latter case, the folder has to also contain .tbi index for each listed file.

`--RegionList`              Local file with regions to test. Has to be in hg38 and in the format 1:100000-12342145. If this is specified, `window` and `Pthresh` arguments are ignored and preselected regions are used in the analysis instead.

`--LiftOver`                Whether input GWAS needs to be lifted over. Values `hg18tohg38` and `hg19tohg38` are allowed. Defaults to `no` (input already in hg38).

**Input format parameters:**

It is possible to specify which column names there are in the input GWAS format. Specify those in `conf/params.conf` or use these flags in the command template.


### Command

Modify this SLURM script with your output paths:

```bash
#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name="RunHyPrColoc"


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
-resume \
-profile singularity,slurm

```

#### Monitoring the job

- After submitting the job, you can monitor it by using:
`watch tail -n 20 slurm-***.out`
- Use `squeue -u [your username]` to see if separate steps of the pipeline run. If only the main job runs for a long time but no steps are in the queue, this might indicate some issue and you should check the slurm log.
- If your job was interrupted due to wall time, just increase the wall time and resubmit it. Pipeline will continue from the step which got interrupted, saving time.
- **When analysis is ready and output chekced, you should remove the directory `work` from your pipeline directory** (this might take considerable disk space), however it means that you need to run all the steps again when resubmitting the job.

### Analysis

After submitting the job, pipeline first constructs the loci based on the GWAS P-value threshold and genomic window. Alternative is to specify tested region(s) in the separate input file. It then runs HyprColoc between all loci for input GWAS and all available/overlapping QTL features in each QTL dataset separately. E.g. if one of your GWAS loci overlaps with 3 genes (gene 1, gene 2, gene 3), HyprColoc is ran between those 4 traits. Output is given for colocalising _trait clusters_, meaning that if your GWAS colocalises with both, gene 1 and gene 2 you also get the posterior probability for GWAS, gene 1, gene 2 colocalisation.

Default settings are used for running HyprColoc command. Pipeline also uses the assumption that all traits are independent which is obviously incorrect for eQTL data (as genes might correlate, GWAS and eQTL might share samples). However, [HyprColoc paper](https://doi.org/10.1038/s41467-020-20885-8) demonstrates that assuming trait independence for correlating traits gives reasonable results.

### Output

Pipeline writes out into output directory two files:

- HyprColocResults.txt: file with colocalisation summaries for all trait clusters (default posterior probability threshold 0.8) for all eQTL files concatenated.
- If `--OutputCsPip TRUE`: 
    
    GwasColocSnpPipResults.txt: file with SNP posterior inclusion probabilities (PIPs) for _colocalisation_ and also if each SNP is part of credible set (95% by default).

Additionally, `pipeline_info` folder is written into output folder, containing Nextflow runtime reports. The most important file from this folder is `Colocalization_report.html`. Download and investigate it, in order to see whether run completed successfully.

### TODO:
- Allow to use other file formats than GWAMA.

### Additional info

You can also run this pipeline on your local computer, without taking advantage of many nodes of HPC. However, it will work for limited analyses. For that you would need:

1. Running instance of Docker on your computer.
2. In your analysis script modify setting `-profile docker`.
3. In the file `conf/base.config` adjust `memory` settings to what is available for your computer and what is allowed to use by Docker. E.g. max 8GB. Hopefully it will work, if you have pre-filtered GWAS file with no more than ~10M SNPs.


## Acknowledgements

This pipeline was written by Urmo Võsa.

HyprColoc is the work of Christopher N. Foley and colleagues:

[Foley, C. N., Staley, J. R., Breen, P. G., Sun, B. B., Kirk, P. D. W., Burgess, S., & Howson, J. M. M. (2021). A fast and efficient colocalization algorithm for identifying shared genetic risk factors across multiple traits. Nature Communications, 12(1), 1–18. https://doi.org/10.1038/s41467-020-20885-8](https://doi.org/10.1038/s41467-020-20885-8)

[Code repo](https://github.com/jrs95/hyprcoloc)

eQTL Catalogue is the work of Nurlan Kerimov and colleagues:

[Kerimov, N., Hayhurst, J. D., Manning, J. R., Walter, P., Kolberg, L., Peikova, K., Samoviča, M., Burdett, T., Jupp, S., Parkinson, H., Papatheodorou, I., Zerbino, D. R., & Alasoo, K. (2020). eQTL Catalogue: a compendium of uniformly processed human gene expression and splicing QTLs. BioRxiv, 2020.01.29.924266. https://doi.org/10.1101/2020.01.29.924266](https://www.biorxiv.org/content/10.1101/2020.01.29.924266v1)

[Website](https://www.ebi.ac.uk/eqtl/)

If you publish results making use of eQTL Catalogue, please cite the paper from Kerimov et al. Please also cite the original eQTL publications and acknowledge the funding they have received, as specificed at the bottom of the page [here](https://www.ebi.ac.uk/eqtl/Studies/).

This tool makes use of Nextflow workflow management tool:
[P. Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017) doi:10.1038/nbt.3820](https://www.nature.com/articles/nbt.3820)