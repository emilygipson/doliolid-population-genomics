# doliolid-population-genomics

Population genomic analyses of *Dolioletta gegenbauri* in the South Atlantic Bight, using low-coverage whole-genome sequencing (lcWGS).

## Analyses

- Mitogenome assembly and characterization
- Mitochondrial diversity and population structure
- Demographic inference
- Nuclear lcWGS
- Comparison of demographic signal between mitochondrial and nuclear data

## Repository structure

```
doliolid-population-genomics/
├── scripts/
│   ├── bash/       # SLURM job scripts for HPC
│   └── R/          # Statistical analysis and visualization
└── markdowns/      # Notes
```

## Software

| Tool | Purpose |
|------|---------|
| ANGSD | Genotype likelihoods, SAF/SFS estimation |
| PCAngsd | PCA from genotype likelihoods |
| MITObim | Iterative mitogenome assembly |
| GetOrganelle | De novo organelle genome assembly |
| BEAST2 | Bayesian Skyline Plot demographic inference |
| BWA, SAMtools | Read mapping and BAM processing |
| R | Statistical analysis and visualization |
