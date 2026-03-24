# doliolid-population-genomics

Population genomic analyses of the doliolid *Dolioletta gegenbauri* in the South Atlantic Bight (SAB), using low-coverage whole-genome sequencing (lcWGS) data.

This repository contains scripts and workflows for investigating genetic structure, reproductive dynamics, and effective population size in a pelagic tunicate characterized by alternating sexual and asexual reproduction and episodic bloom formation.

## Overview

- **SNP discovery and genotype likelihood estimation** from lcWGS data using ANGSD
- **Population structure analysis** via PCA (PCAngsd) and admixture
- **Imputation and phasing** with BEAGLE
- **Forward simulations** (SLiM) modeling bloom dynamics and mixed reproductive strategies
- **Approximate Bayesian Computation (ABC)** to evaluate demographic scenarios and estimate selfing rates
- **Reference transcriptome assembly** for *D. gegenbauri* using RNA-seq data

## Software and dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| ANGSD | | Genotype likelihoods, SAF/SFS estimation, SNP calling |
| PCAngsd | | PCA, admixture, selection scans from genotype likelihoods |
| BEAGLE | | Genotype imputation and phasing |
| SLiM | | Forward genetic simulations |
| Trinity | | De novo transcriptome assembly |
| R | >= 4.x | Statistical analysis, visualization, ABC |

## Data availability

Raw sequence data are archived at NCBI SRA (accession pending). Large data files are not included in this repository.

## Associated publications

Gipson, E.E. et al. *In preparation.*

## Contact

Emily Gipson - Emily.Gipson@uga.edu | GitHub: @emilygipson
