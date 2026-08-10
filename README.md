# Rhizosphere Community Modeling

A reproducible bioinformatics pipeline for reconstructing metagenome-assembled genomes (MAGs), generating genome-scale metabolic models (GEMs), and performing community metabolic simulations from rhizosphere metagenomic data.

---

## Overview

This repository contains the scripts used to process rhizosphere metagenomic sequencing data, reconstruct microbial genomes, generate genome-scale metabolic models, simulate microbial communities using MICOM, and perform downstream statistical analyses.

The workflow was developed as part of an M.Sc. thesis at the University of Potsdam.

---

## Workflow

The pipeline consists of four major stages:

### 1. MAG reconstruction

* Download sequencing data
* Quality control
* Adapter and quality trimming
* Host read removal
* Metagenome assembly
* Contig filtering
* Read mapping
* Coverage estimation
* Metagenomic binning
* Bin refinement
* Quality assessment
* High-quality MAG extraction
* Taxonomic classification
* Dereplication
* Gene prediction

Scripts:

```text
scripts/MAGs/
```

---

### 2. Genome-scale metabolic modeling

* CarveMe model reconstruction
* MICOM community simulations
* Flux variability analysis (FVA)
* Memote quality assessment
* GEM quality control

Scripts:

```text
scripts/GEMs/
```

---

### 3. Statistical analyses

Downstream analyses include:

* Differential flux analysis
* Flux–trait correlations
* Principal component analysis (PCA)
* Treatment-adjusted regression analyses
* Plant trait comparisons
* Publication-quality figures

Scripts:

```text
scripts/Statistics/
```

---

## Repository structure

```text
rhizosphere-community-modeling/
├── README.md
├── environments/
├── example_data/
└── scripts/
    ├── MAGs/
    ├── GEMs/
    └── Statistics/
```

---

## Software environments

Separate Conda environments are provided for different stages of the workflow.

| Environment    | Purpose                                            |
| -------------- | -------------------------------------------------- |
| `megahit_env`  | Metagenome assembly                                |
| `binning-env`  | Metagenomic binning and refinement                 |
| `gtdbtk-2.1.1` | Taxonomic classification                           |
| `prokka_env`   | Genome annotation                                  |
| `carveme_env`  | Genome-scale metabolic model reconstruction        |
| `micom_env`    | Community metabolic modeling                       |
| `corr_env`     | Statistical analyses and visualization             |
| `memote_env`   | GEM quality assessment                             |
| `qc_gems`      | Additional metabolic model quality-control scripts |
| `seqkit_env`   | FASTA/FASTQ processing                             |

Environment files are located in:

```text
environments/
```

To create an environment:

```bash
conda env create -f environments/micom_env.yml
```

Activate it using:

```bash
conda activate micom_env
```

---

## HPC requirements

The pipeline was developed and tested on a SLURM-managed high-performance computing (HPC) cluster running Linux.

Several scripts require software provided through environment modules, including tools such as:

* BWA
* SAMtools
* BioPerl

Adjust module names as required for your computing environment.

---

## Input data

The pipeline expects paired-end Illumina metagenomic reads as input.

Example input files are provided in:

```text
example_data/
```

Large sequencing datasets are intentionally excluded from this repository.

---

## Output

The pipeline generates:

* Metagenome assemblies
* MAGs
* Dereplicated genomes
* Genome-scale metabolic models
* MICOM community models
* Exchange fluxes
* Flux variability analysis results
* Statistical summaries
* Publication-quality figures

---

## Citation

If you use this workflow, please cite the corresponding thesis and the software packages used throughout the analysis, including MEGAHIT, MetaBAT2, MaxBin2, CONCOCT, DAS Tool, GTDB-Tk, Prokka, CarveMe, MICOM, CheckM, dRep, and Memote.

---

## License

This project is licensed under the MIT License. See the LICENSE file for details.

---

## Contact

Mihriban Seyis

University of Potsdam
