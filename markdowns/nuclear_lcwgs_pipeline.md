# Nuclear lcWGS Pipeline

Low-coverage whole genome sequencing pipeline for 118 *Dolioletta gegenbauri* individuals. Reads are mapped to a de novo transcriptome reference and analyzed entirely within the ANGSD genotype likelihood framework — no hard genotype calls anywhere in the pipeline.

---

## Overview

We have paired-end Illumina reads from 120 individuals. Two (DL5_02 and DL5_04) were dropped for extremely low mapped read counts (3,803 and 12,145), leaving 118 for all nuclear analyses. Coverage is ~0.2x per sample, so everything downstream uses genotype likelihoods.

The reference is a de novo transcriptome (88,916 CDS contigs, 93.6% BUSCO completeness), not a genome assembly, so all analyses are restricted to coding regions and contigs are short.

---

## Step 1: Read mapping and duplicate removal

**Scripts:**
- [`scripts/bash/20260319_map_markdup_21.sh`](../scripts/bash/20260319_map_markdup_21.sh) — BWA-MEM mapping + SAMtools sort + Picard MarkDuplicates
- [`scripts/bash/20260319_addRG_markdup_21.sh`](../scripts/bash/20260319_addRG_markdup_21.sh) — adds read groups and re-runs MarkDuplicates

Paired-end reads are mapped with BWA-MEM (default settings), coordinate-sorted, and deduplicated with Picard MarkDuplicates (`REMOVE_DUPLICATES=true` — physically removes duplicates, doesn't just flag them). Read groups are added with AddOrReplaceReadGroups. These scripts handle the final 21 samples added to bring the dataset to 120; the earlier 99 were processed with an equivalent pipeline.

Output is one deduplicated, indexed BAM per sample.

---

## Step 2: Coverage assessment and downsampling

**Scripts:**
- [`scripts/bash/20260319_downsample_check.sh`](../scripts/bash/20260319_downsample_check.sh) — mapped read counts per sample
- [`scripts/bash/20260319_build_bamlist_depth.sh`](../scripts/bash/20260319_build_bamlist_depth.sh) — builds the BAM list for ANGSD

Mapped read counts come from `samtools idxstats`. Dataset median is 186,635 mapped reads. Fourteen samples with much higher coverage were downsampled to the median with `samtools view -s` (seeded for reproducibility) to avoid coverage-driven artifacts in PCA and diversity estimates.

The BAM list (one path per line, consistent ordering) is the input for everything downstream.

---

## Step 3: SAF likelihoods and SFS

**Script:**
- [`scripts/bash/20260320_angsd_saf_sfs_118.sh`](../scripts/bash/20260320_angsd_saf_sfs_118.sh)

ANGSD estimates site allele frequency (SAF) likelihoods across the reference, then `realSFS` computes the global folded SFS. Diversity and neutrality stats come from the SFS via `thetaStat`.

**Filters:**
- `-GL 1` (SAMtools GL model)
- `-minMapQ 30`, `-minQ 20`, `-baq 1`
- `-setMinDepth 2`, `-setMaxDepth 50` — the max depth filter is important here because paralogous transcripts mapping to the same contig inflate coverage and fake heterozygosity
- `-minInd 24` (~20% of 118) — can't be too strict at 0.2x or you lose most sites
- `-fold 1` (no outgroup for polarization)

**Output:** pi = 0.004462, theta_W = 0.004750, Tajima's D ~ -0.30 (not significant).

---

## Step 4: Genotype likelihoods for PCA

**Script:**
- [`scripts/bash/20260320_angsd_beagle_118.sh`](../scripts/bash/20260320_angsd_beagle_118.sh)

SNP calling and GL output in beagle format (three columns per individual per SNP: P(RR), P(RA), P(AA)). Same quality filters as Step 3, plus:

- `-SNP_pval 1e-6`
- `-minMaf 0.05` — rare variants add noise to PCA
- `-doGlf 2` (beagle output)
- `-doMajorMinor 1`, `-doMaf 1`

67,091 SNPs pass filters.

---

## Step 5: PCA

**Script:**
- [`scripts/bash/20260321_pcangsd_118_v2.sh`](../scripts/bash/20260321_pcangsd_118_v2.sh)

PCAngsd v1.2 on the beagle file, `--iter 500`, `--maf_iter 500`.

MAP test finds one significant PC (PC1, 5.7% variance). PC2–PC10 are each <0.9% variance.

Eigenvalues and eigenvectors are plotted in R, colored by bloom/non-bloom, year, and life stage.

---

## Step 6: Bloom vs. non-bloom Fst

**Script:**
- [`scripts/bash/20260320_bloom_fst_118.sh`](../scripts/bash/20260320_bloom_fst_118.sh)

Separate SAF estimation for bloom (n = 80) and non-bloom (n = 38) using the same filters, then 2D folded SFS estimation with `realSFS`, then Fst.

Weighted Fst = 0.002, unweighted = 0.007.

---

## Software versions

| Tool | Version |
|------|---------|
| BWA-MEM | 0.7.18 |
| SAMtools | 1.21 |
| Picard | 3.3.0 |
| ANGSD | 0.940 |
| PCAngsd | 1.2 |

---

## Notes

All scripts are SLURM jobs for the UGA GACRC Sapelo2 cluster. Module loads and SBATCH headers are cluster-specific but the tool invocations are portable.

NgsRelate (`20260321_ngsrelate_118_v2.sh`) was attempted but returned zero overlapping sites for all 6,903 pairwise comparisons — a limitation of extremely low coverage on short CDS contigs.
