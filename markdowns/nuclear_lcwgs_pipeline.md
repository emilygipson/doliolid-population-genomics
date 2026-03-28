# Nuclear lcWGS Pipeline

Low-coverage whole genome sequencing pipeline for 118 *Dolioletta gegenbauri* individuals mapped to a de novo transcriptome reference. All analyses use genotype likelihoods (ANGSD framework), no hard genotype calls.

---

## Overview

120 individuals, 2 dropped (DL5_02 and DL5_04, low mapped read counts), 118 retained. ~0.2x coverage per sample. Reference: de novo transcriptome, 88,916 CDS contigs, 93.6% BUSCO completeness.

---

## Step 1: Read mapping and duplicate removal

**Scripts:**
- [`scripts/bash/20260319_map_markdup_21.sh`](../scripts/bash/20260319_map_markdup_21.sh)
- [`scripts/bash/20260319_addRG_markdup_21.sh`](../scripts/bash/20260319_addRG_markdup_21.sh)

BWA-MEM (default settings) → SAMtools sort → Picard MarkDuplicates (`REMOVE_DUPLICATES=true`) → AddOrReplaceReadGroups. These scripts handle the final 21 samples; the earlier 99 were processed with an equivalent pipeline.

Output: one deduplicated, indexed BAM per sample.

---

## Step 2: Coverage assessment and downsampling

**Scripts:**
- [`scripts/bash/20260319_downsample_check.sh`](../scripts/bash/20260319_downsample_check.sh)
- [`scripts/bash/20260319_build_bamlist_depth.sh`](../scripts/bash/20260319_build_bamlist_depth.sh)

Mapped read counts from `samtools idxstats`. Median: 186,635. Fourteen high-coverage samples downsampled to median with `samtools view -s` (seeded). BAM list built for downstream ANGSD input.

---

## Step 3: SAF likelihoods and SFS

**Script:**
- [`scripts/bash/20260320_angsd_saf_sfs_118.sh`](../scripts/bash/20260320_angsd_saf_sfs_118.sh)

SAF likelihoods → `realSFS` → folded SFS → `thetaStat` for pi and theta_W.

**Filters:** `-GL 1`, `-minMapQ 30`, `-minQ 20`, `-baq 1`, `-setMinDepth 2`, `-setMaxDepth 50`, `-minInd 24`, `-fold 1`

**Output:** pi = 0.004462, theta_W = 0.004750, Tajima's D ~ -0.30.

---

## Step 4: Genotype likelihoods for PCA

**Script:**
- [`scripts/bash/20260320_angsd_beagle_118.sh`](../scripts/bash/20260320_angsd_beagle_118.sh)

Same quality filters as Step 3, plus: `-SNP_pval 1e-6`, `-minMaf 0.05`, `-doGlf 2`, `-doMajorMinor 1`, `-doMaf 1`

67,091 SNPs pass filters.

---

## Step 5: PCA

**Script:**
- [`scripts/bash/20260321_pcangsd_118_v2.sh`](../scripts/bash/20260321_pcangsd_118_v2.sh)

PCAngsd v1.2, `--iter 500`, `--maf_iter 500`. MAP test: 1 significant PC (PC1, 5.7% variance). PC2–PC10 each <0.9%.

Eigenvalues and eigenvectors plotted in R.

---

## Step 6: Bloom vs. non-bloom Fst

**Script:**
- [`scripts/bash/20260320_bloom_fst_118.sh`](../scripts/bash/20260320_bloom_fst_118.sh)

Separate SAF for bloom (n = 80) and non-bloom (n = 38) → 2D folded SFS → `realSFS fst`.

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

NgsRelate (`20260321_ngsrelate_118_v2.sh`) was attempted but returned zero overlapping sites for all 6,903 pairwise comparisons.
