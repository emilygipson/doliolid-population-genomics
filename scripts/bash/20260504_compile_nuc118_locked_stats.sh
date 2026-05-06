#!/bin/bash
#SBATCH --job-name=compile_nuc118_locked_stats
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu
#SBATCH --output=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/compile_nuc118_locked_stats_%j.out
#SBATCH --error=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/compile_nuc118_locked_stats_%j.err

# Compile all locked nuclear (n=118) summary statistics into a single CSV
# with file-level provenance. Every value is recomputed from raw output
# files; nothing is hardcoded from documentation.
#
# Run any time to verify locked values have not shifted, or to regenerate
# the CSV after a rerun. Self-contained: no external interpreter required.

set -euo pipefail

module load angsd/0.940-GCC-12.3.0
module load R/4.4.1-gfbf-2023b

WORKDIR=/scratch/eeg37520/nuc_120_mitomatch_redo
PESTPG=${WORKDIR}/tajimaD_20260502/global_perContig.pestPG
GLOBAL_SFS=${WORKDIR}/angsd/nuc118_global.sfs
BEAGLE=${WORKDIR}/angsd/nuc118_beagle_minInd24.beagle.gz
PCANGSD_COV=${WORKDIR}/pcangsd/nuc118_pcangsd_minInd24_v2.cov

OUTDIR=${WORKDIR}/locked_stats
mkdir -p ${OUTDIR}
OUT=${OUTDIR}/nuc118_locked_summary_stats_$(date +%Y%m%d).csv

DATE_RECOMPUTED=$(date +%Y-%m-%d)

# CSV header
echo "metric,value,n_samples,total_sites,source_file,formula,date_locked,date_recomputed" > ${OUT}

echo "=== File presence check ==="
for F in ${PESTPG} ${GLOBAL_SFS} ${BEAGLE} ${PCANGSD_COV}; do
    if [ -f "${F}" ]; then
        echo "FOUND:    ${F}"
    else
        echo "MISSING:  ${F}"
    fi
done
echo

# --- Global Tajima's D, pi, theta_W from pestPG ---
if [ -f "${PESTPG}" ]; then
    D=$(awk 'NR>1 {sum_d += $9 * $14; sum_n += $14} END {printf "%.6f", sum_d/sum_n}' ${PESTPG})
    NSITES=$(awk 'NR>1 {sum_n += $14} END {print sum_n}' ${PESTPG})
    NCONTIGS=$(awk 'NR>1 {n++} END {print n}' ${PESTPG})
    PI=$(awk 'NR>1 {sum_tp += $5; sum_n += $14} END {printf "%.6f", sum_tp/sum_n}' ${PESTPG})
    THETA_W=$(awk 'NR>1 {sum_tw += $4; sum_n += $14} END {printf "%.6f", sum_tw/sum_n}' ${PESTPG})

    echo "tajima_D_global,${D},118,${NSITES},${PESTPG},sum(D_i*nSites_i)/sum(nSites_i) per-contig pestPG cols 9 and 14,2026-05-02,${DATE_RECOMPUTED}" >> ${OUT}
    echo "pi_global,${PI},118,${NSITES},${PESTPG},sum(tP)/sum(nSites) per-site pi pestPG cols 5 and 14,2026-03-21,${DATE_RECOMPUTED}" >> ${OUT}
    echo "theta_W_global,${THETA_W},118,${NSITES},${PESTPG},sum(tW)/sum(nSites) per-site Watterson theta pestPG cols 4 and 14,2026-03-21,${DATE_RECOMPUTED}" >> ${OUT}
    echo "n_contigs,${NCONTIGS},118,${NSITES},${PESTPG},rows in pestPG excluding header,2026-03-21,${DATE_RECOMPUTED}" >> ${OUT}
    echo "n_sites_analyzed,${NSITES},118,${NSITES},${PESTPG},sum of nSites column,2026-03-21,${DATE_RECOMPUTED}" >> ${OUT}
fi

# --- Global folded SFS file marker ---
if [ -f "${GLOBAL_SFS}" ]; then
    SFS_LEN=$(awk '{print NF}' ${GLOBAL_SFS})
    echo "sfs_global_length,${SFS_LEN},118,NA,${GLOBAL_SFS},folded SFS entry count expected 237 for 118 diploids,2026-03-21,${DATE_RECOMPUTED}" >> ${OUT}
fi

# --- SNP count for PCA ---
if [ -f "${BEAGLE}" ]; then
    NSNP=$(zcat ${BEAGLE} | tail -n +2 | wc -l)
    echo "n_snps_pca,${NSNP},118,NA,${BEAGLE},rows in beagle.gz minus header SNP_pval 1e-6 minMaf 0.05,2026-03-21,${DATE_RECOMPUTED}" >> ${OUT}
fi

# --- PCA variance explained from PCAngsd .cov via R ---
if [ -f "${PCANGSD_COV}" ]; then
    PC_VARS=$(Rscript --vanilla -e "
.libPaths(c('/scratch/eeg37520/Rlibs', .libPaths()))
cov <- as.matrix(read.table('${PCANGSD_COV}'))
e <- eigen(cov)
pcts <- 100 * e\$values / sum(e\$values)
cat(sprintf('%.4f', pcts[1]), sprintf('%.4f', pcts[2]), sep=',')
" 2>/dev/null || echo "NA,NA")
    PC1=$(echo ${PC_VARS} | cut -d',' -f1)
    PC2=$(echo ${PC_VARS} | cut -d',' -f2)
    echo "pca_pc1_variance_pct,${PC1},118,NA,${PCANGSD_COV},eigenvalue[1]/sum(eigenvalues)*100,2026-03-21,${DATE_RECOMPUTED}" >> ${OUT}
    echo "pca_pc2_variance_pct,${PC2},118,NA,${PCANGSD_COV},eigenvalue[2]/sum(eigenvalues)*100,2026-03-21,${DATE_RECOMPUTED}" >> ${OUT}
fi

# --- Bloom-vs-non-bloom Fst (cruise-level, March 21; SUPERSEDED) ---
OLD_FST_IDX=${WORKDIR}/fst/bloom_nonbloom_fst.fst.idx
if [ -f "${OLD_FST_IDX}" ]; then
    OLD_FST_TXT=${OUTDIR}/old_cruise_fst_stats_recomputed.txt
    realSFS fst stats ${OLD_FST_IDX} > ${OLD_FST_TXT} 2>&1 || true
    OLD_FST_LINE=$(tail -1 ${OLD_FST_TXT})
    OLD_FST_W=$(echo "${OLD_FST_LINE}" | awk '{print $NF}')
    OLD_FST_U=$(echo "${OLD_FST_LINE}" | awk '{print $(NF-1)}')
    echo "fst_bloom_v_nonbloom_unweighted_CRUISE_LEVEL_SUPERSEDED,${OLD_FST_U},118 (80/38),NA,${OLD_FST_IDX},realSFS fst stats CRUISE-LEVEL PARTITION superseded,2026-03-21,${DATE_RECOMPUTED}" >> ${OUT}
    echo "fst_bloom_v_nonbloom_weighted_CRUISE_LEVEL_SUPERSEDED,${OLD_FST_W},118 (80/38),NA,${OLD_FST_IDX},realSFS fst stats CRUISE-LEVEL PARTITION superseded,2026-03-21,${DATE_RECOMPUTED}" >> ${OUT}
fi

# --- Bloom-vs-non-bloom Fst (station-level, May 4 rerun) ---
NEW_FST_DIR=${WORKDIR}/fst_stationlevel_20260504
for FLAG in 0 1; do
    STATS=${NEW_FST_DIR}/bloom_nonbloom_stationlevel_fst_whichFst${FLAG}_stats.txt
    if [ -f "${STATS}" ]; then
        FST_LINE=$(tail -1 ${STATS})
        FST_W=$(echo "${FST_LINE}" | awk '{print $NF}')
        FST_U=$(echo "${FST_LINE}" | awk '{print $(NF-1)}')
        echo "fst_bloom_v_nonbloom_unweighted_stationlevel_whichFst${FLAG},${FST_U},118 (63/55),NA,${STATS},realSFS fst stats station-level partition whichFst ${FLAG},2026-05-04,${DATE_RECOMPUTED}" >> ${OUT}
        echo "fst_bloom_v_nonbloom_weighted_stationlevel_whichFst${FLAG},${FST_W},118 (63/55),NA,${STATS},realSFS fst stats station-level partition whichFst ${FLAG},2026-05-04,${DATE_RECOMPUTED}" >> ${OUT}
    fi
done

# --- Print summary ---
echo
echo "=== Locked nuclear n=118 stats ==="
column -t -s, ${OUT}
echo
echo "Saved: ${OUT}"
