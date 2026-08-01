# Nuclear lcWGS workflow

ANGSD steps for the *D. gegenbauri* nuclear low-coverage data, n=118. SAF, folded
SFS, thetas, pairwise Fst, PCAngsd. Mapping and dedup happen upstream in
scripts/bash/ and aren't part of this.

Full pipeline writeup: markdowns/nuclear_lcwgs_pipeline.md

## Groups and pairs

Groups are bamlists, named in config.yaml. Every per-group rule runs across all
of them. Pairs name two groups and get Fst under both -whichFst 0 (Reynolds)
and -whichFst 1 (Bhatia).

    bamlist -> saf.idx -> sfs -> thetas.idx -> pestPG
    bamlist -> saf.idx -> 2d.sfs -> fst.idx -> stats.txt
    bamlist -> beagle.gz -> cov

To add a partition: drop in a bamlist, add one line to config.yaml.

## Running

    conda deactivate
    module load snakemake/8.27.0-foss-2024a
    snakemake --profile profiles/slurm -n
    snakemake --profile profiles/slurm

Deactivate conda first. With base active it still works, but only because of
path ordering, which isn't worth relying on.

The regions entry in config.yaml limits ANGSD to a contig subset. Useful for
testing. Blank it for the full 88,916.

## Parameters

-setMaxDepth 300, -setMinDepth 2, -minMapQ 30, -minQ 20, -baq 1. -minInd is
20% of group size, except the 118-sample set, which is 24 (20% truncates to 23).

Reference is the one-CDS-per-gene transcriptome at
/scratch/eeg37520/transdecoder_homology/cdhit/one_per_gene/. Needs a .fai.

## Caveats

The bloom/non-bloom bamlists here are the May 2026 station-level split. That
partition is being redone. Fst computed under different partitions isn't
comparable.

SAF on 118 samples genome-wide is the slow step. Everything after it is fast.
