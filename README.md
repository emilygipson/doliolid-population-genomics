# doliolid-population-genomics

Population genomic analyses of the doliolid tunicate *Dolioletta gegenbauri* in the South Atlantic Bight, using low-coverage whole-genome sequencing (lcWGS) 

## Overview

This repository contains scripts and analysis workflows for investigating genetic structure, bloom dynamics, and mito-nuclear discordance in a bloom-forming pelagic tunicate.

Key analyses include:

- **Mitogenome assembly and characterization** 
- **Population diversity** — haplotype diversity, nucleotide diversity, site frequency spectrum, neutrality tests, per-gene dN/dS
- **Bloom vs non-bloom comparisons** — Hudson's Snn, AMOVA, Fst, clone detection
- **Bayesian Skyline Plot** demographic inference (BEAST2) with mutation rate sensitivity analysis
- **Nuclear lcWGS analysis** — genotype likelihoods via ANGSD, PCA via PCAngsd, SFS-based diversity and Fst, BEAGLE imputation
- **Mito-nuclear discordance** — comparison of demographic signals across genomes

## Repository structure
```
doliolid-population-genomics/
├── scripts/
│   ├── bash/       # SLURM job scripts for HPC (Sapelo2)
│   └── R/          # Statistical analysis and visualization
├── markdowns/      # Step-by-step analysis documentation
└── figures/        # Selected output figures
```

## Software and dependencies

| Tool | Purpose |
|------|---------|
| ANGSD | Genotype likelihoods, SAF/SFS estimation |
| PCAngsd | PCA and admixture from genotype likelihoods |
| BEAGLE | Genotype imputation and phasing |
| MITObim | Iterative mitogenome assembly |
| GetOrganelle | De novo organelle genome assembly |
| BEAST2 | Bayesian Skyline Plot demographic inference |
| FreeBayes | Variant calling for mitogenome validation |
| SAMtools / BWA | Read mapping and BAM processing |
| R | Statistical analysis, visualization, AMOVA, clone detection |



## Contact

Emily Gipson — Emily.Gipson@uga.edu
GitHub: [@emilygipson](https://github.com/emilygipson)
