# Genotyping by sequencing (DNA- or RNA-seq data)

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.0-brightgreen.svg)](https://snakemake.bitbucket.io)    
[![Miniconda](https://img.shields.io/badge/miniconda-blue.svg)](https://conda.io/miniconda)
[![GitHub
Issues](https://img.shields.io/github/issues/bleekerlab/snakemake_rnaseq.svg)](https://github.com/bleekerlab/snakemake_rnaseq/issues)
[![Documentation Status](https://readthedocs.org/projects/snakemake-rnaseq-pipeline/badge/?version=latest)](https://snakemake-rnaseq-pipeline.readthedocs.io/en/latest/?badge=latest)

<!-- MarkdownTOC autolink="true" levels="1,2" -->

- [Description](#description)
	- [Description](#description-1)
	- [Input files](#input-files)
	- [Output files](#output-files)
	- [Prerequisites: what you should know before using this pipeline](#prerequisites-what-you-should-know-before-using-this-pipeline)
	- [Content of this GitHub repository](#content-of-this-github-repository)
- [Installation and usage \(local machine\)](#installation-and-usage-local-machine)
	- [Installation](#installation)
	- [Usage](#usage)
	- [Configuration :pencil2:](#configuration-pencil2)
	- [Dry run](#dry-run)
	- [Real run](#real-run)
- [Installation and usage \(HPC cluster\)](#installation-and-usage-hpc-cluster)
	- [Installation](#installation-1)
	- [SLURM sbatch usage](#slurm-sbatch-usage)
	- [SLURM interactive usage](#slurm-interactive-usage)
	- [Retrieving your results](#retrieving-your-results)
	- [Useful links](#useful-links)
- [Directed Acyclic Graph of jobs](#directed-acyclic-graph-of-jobs)
- [References :green_book:](#references-green_book)
	- [Authors](#authors)
	- [Pipeline dependencies](#pipeline-dependencies)
	- [Acknowledgments :clap:](#acknowledgments-clap)
- [Citation](#citation)

<!-- /MarkdownTOC -->


# Description

A Snakemake pipeline that calls SNPs (Single Nucleotide Polymorphisms) from either DNA-seq or mRNA-seq data. It trims and aligns DNA/mRNA-seq fastq to a reference genome, then call SNPs before computing the number of SNPs per genomic window (e.g. per 1Mb). 
This pipeline can process single or paired-end data and is mostly suited for Illumina sequencing data. 

## Steps

1. The raw fastq files will be trimmed for adaptors and quality checked with `fastp`.  
2. The genome sequence FASTA file will be used for the mapping step of the trimmed reads using either `bwa` or `STAR`. 
3. The reference genome will be sliced into bins of a pre-determined size (e.g. 1Mb).  
4. SNPs are called based on the alignment `.bam` files generating one VCF file per sample.
5. SNPs are filtered based on a read depth threshold and SNP quality. 
6. The obtained filtered VCF file(s) are converted to the BED format `.bed` to allow computation. 
7. The number of overlapping SNPs per genome bin is computed using BEDOPS `bedmap` operator. 


## Input files
* __DNA or RNA-seq fastq files__ as listed in the `config/samples.tsv` file. Specify a sample name (e.g. "Sample_A") in the `sample` column and the paths to the forward read (`fq1`) and to the reverse read (`fq2`). If you have single-end reads, leave the `fq2` column empty. 
* __A genomic reference in FASTA format__. For instance, a fasta file containing the 12 chromosomes of tomato (*Solanum lycopersicum*).
* __(for RNA-seq datasets) A genome annotation file in the [GTF format](https://useast.ensembl.org/info/website/upload/gff.html)__. You can convert a GFF annotation file format into GTF with the [gffread program from Cufflinks](http://ccb.jhu.edu/software/stringtie/gff.shtml): `gffread my.gff3 -T -o my.gtf`. :warning: for featureCounts to work, the _feature_ in the GTF file should be `exon` while the _meta-feature_ has to be `transcript_id`. This will be converted to a BED file to compute the number of genes per genomic bin. 

Below is an example of a GTF file format. :warning: a real GTF file does not have column names (seqname, source, etc.). Remove all non-data rows. 

| seqname | source | feature | start  | end  | score | strand | frame | attributes |
|-----------|------------|------|------|------|---|---|---|----------------------------------------------------------------------------------------------------|
| SL4.0ch01 | maker_ITAG | CDS  | 279  | 743  | . | + | 0 | transcript_id "Solyc01g004000.1.1"; gene_id "gene:Solyc01g004000.1"; gene_name "Solyc01g004000.1"; |
| SL4.0ch01 | maker_ITAG | exon | 1173 | 1616 | . | + | . | transcript_id "Solyc01g004002.1.1"; gene_id "gene:Solyc01g004002.1"; gene_name "Solyc01g004002.1"; |
| SL4.0ch01 | maker_ITAG | exon | 3793 | 3971 | . | + | . | transcript_id "Solyc01g004002.1.1"; gene_id "gene:Solyc01g004002.1"; gene_name "Solyc01g004002.1"; |

## Output files

* __One .tsv file per sample summarising the number of counts per genome bin__ called `[sample name].counts.tsv`. This table is used for plotting. 
* __fastp QC reports__: one per fastq file.

## Prerequisites: what you should know before using this pipeline
- Some command of the Unix Shell to connect to a remote server where you will execute the pipeline. You can find a good tutorial from the Software Carpentry Foundation [here](https://swcarpentry.github.io/shell-novice/) and another one from Berlin Bioinformatics [here](http://bioinformatics.mdc-berlin.de/intro2UnixandSGE/unix_for_beginners/README.html).
- Some command of the Unix Shell to transfer datasets to and from a remote server (to transfer sequencing files and retrieve the results/). The Berlin Bioinformatics Unix begginer guide available [here](http://bioinformatics.mdc-berlin.de/intro2UnixandSGE/unix_for_beginners/README.html)) should be sufficient for that (check the `wget` and `scp` commands).
- An understanding of the steps of a canonical RNA-Seq analysis (trimming, alignment, etc.). You can find some info [here](https://bitesizebio.com/13542/what-everyone-should-know-about-rna-seq/).

## Content of this GitHub repository
- `Snakefile`: a master file that contains the desired outputs and the rules to generate them from the input files.
- `config/samples.tsv`:  a file containing sample names and the paths to the forward and eventually reverse reads (if paired-end). **This file has to be adapted to your sample names before running the pipeline**.
- `config/config.yaml`: the configuration files making the Snakefile adaptable to any input files, genome and parameter for the rules.
- `config/refs/`: a folder containing
  - a genomic reference in fasta format. The `S_lycopersicum_chromosomes.4.00.chrom1.fa` is placed for testing purposes.
  - a GTF annotation file. The `ITAG4.0_gene_models.sub.gtf` for testing purposes.
- The `environment.yaml` is used by the conda package manager to create a working environment (see below).

You have to create the `config/fastq/` folder that should contain test fastq files. These test fastq files are subsetted paired-end fastq files used to test locally the pipeline. They are generated using [Seqtk](https://github.com/lh3/seqtk):`seqtk sample -s100 <inputfile> 250000 > <output file>`. 
Files can be downloaded from [Zenodo](https://doi.org/10.5281/zenodo.4085315).  

# Installation and usage (local machine)

## Installation

You will need a local copy of the GitHub `genotype_by_sequencing` repository on your machine. You can either:
- use git in the shell: `git clone git@github.com:BleekerLab/genotype_by_sequencing.git`.
- click on ["Clone or download"](https://github.com/BleekerLab/genotype_by_sequencing/archive/master.zip) and select `download`. 
- Then navigate inside the `genotype_by_sequencing` folder using Shell commands.

## Usage 

## Configuration :pencil2:
You'll need to change a few things to accomodate this pipeline to your needs. 

1. Make sure you have changed the parameters in the `config/config.yaml` file that specifies where to find the sample data file, the genomic and transcriptomic reference fasta files to use and the parameters for certains rules etc. __in particular, indicate whether your data are DNA-seq or RNA-seq data__  
The type of input data has to be either `DNA` or `RNA` in the `datatype` section of the `config.yaml`. 

This file is used so the `Snakefile` does not need to be changed when locations or parameters need to be changed.

### :round_pushpin: conda
Using the conda package manager, you need to create an environment where core softwares such as `Snakemake` will be installed.   
1. Install the [Miniconda3 distribution (>= Python 3.7 version)](https://docs.conda.io/en/latest/miniconda.html) for your OS (Windows, Linux or Mac OS X).  
2. Inside a Shell window (command line interface), create a virtual environment named `gbs` using the `environment.yaml` file with the following command: `conda install -c conda-forge mamba --yes && mamba env create --file environment.yaml`. The `mamba` package manager is faster than conda. 
3. Then, before you run the Snakemake pipeline, activate this virtual environment with `conda activate gbs`.

## Dry run
- Use the `snakemake -np` to perform a dry run that prints out the rules and commands.

## Real run
With conda: `snakemake --cores 10`


# Installation and usage (HPC cluster)

## Installation
You will need a local copy of the GitHub `genotype_by_sequencing` repository on your machine. On a HPC system, you will have to clone it using the Shell command-line: `git clone git@github.com:BleekerLab/genotype_by_sequencing.git`.
- click on ["Clone or download"](https://github.com/BleekerLab/genotype_by_sequencing/archive/master.zip) and select `download`. 
- Then navigate inside the `genotype_by_sequencing` folder using Shell commands.

## SLURM sbatch usage

See the detailed protocol [here](./hpc/README.md). 

Here is an example script to be saved as `my_run.sh` and executed with SLURM as `sbatch my_run.sh`  
This will submit the batch script to SLURM sbatch. 

```bash
#!/bin/bash
#
#SBATCH --job-name=snakemake_gbs        # job name
#SBATCH  --time=24:00:00                # mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem-per-cpu=8G                # RAM requested per job (in Snakemake, one rule = one job)
#SBATCH --output=parallel_%j.log        # standard output and error log (%j substitutes the JOB ID)


#SBATCH --nodes=1                       # Run all processes on a single node	
#SBATCH --cpus-per-task=30              # Number of CPUs per task

source activate genotype

srun snakemake -j $SLURM_CPUS_PER_TASK
```

## SLURM interactive usage 

To start an interactive session, do: 

- Step1: `srun --time=24:00:00 --mem-per-cpu=8G --cpus-per-task=10 --pty bash -i`    
- Step2: `conda activate gbs`    
- Step3: `snakemake -j 10` (since we specified 10 cpus per task)    
- Step4: `exit`  

This starts an interactive bash with 10 CPUs allocated per task and 8G of RAM per CPU. 

## Retrieving your results

For instance, download the results but exclude bam files (too big).   
`rsync -a -v -e ssh --exclude="*bam"  mgallan1@omics-h0.science.uva.nl:/zfs/omics/personal/mgallan1/workspace/genotype_by_sequencing/results/ [local directory]`

## Useful links
- [Snakemake with SLURM](https://accio.github.io/programming/2020/06/16/Snakemake-with-slurm.html)
- [`sbatch` manual with its options](https://slurm.schedmd.com/sbatch.html)


# Directed Acyclic Graph of jobs
![dag](./dag.png)

# References :green_book:

## Authors
- Marc Galland, m.galland@uva.nl 
- Tijs Bliek, m.bliek@uva.nl


## Pipeline dependencies
* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [fastp](https://github.com/OpenGene/fastp)
* [STAR](https://github.com/alexdobin/STAR)   
* [bcftools](https://samtools.github.io/bcftools/howtos/index.html)  
* [bedtools](https://bedtools.readthedocs.io/en/latest/index.html)
* [bedops](https://bedops.readthedocs.io/en/latest/index.html)


## Acknowledgments :clap:
[Johannes Köster](https://johanneskoester.bitbucket.io/); creator of Snakemake. 

# Citation
If you use this software, please use the following citation:  

Bliek T. and Galland M. (2021). Genotyping by sequencing pipeline (version 0.1.0). 