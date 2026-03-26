#!/bin/bash
#SBATCH --job-name=dgeg_freebayes
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/eeg37520/doliolid_popgen/mitobim/vcf_dgeg/logs/freebayes_%j.out
#SBATCH --error=/scratch/eeg37520/doliolid_popgen/mitobim/vcf_dgeg/logs/freebayes_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu

# Joint variant calling on 105 samples remapped to D. gegenbauri reference
# Run AFTER remap array job completes
# Date: 2026-03-10
# Fix: corrected freebayes module, absolute log paths, added email

module load freebayes/1.3.7-gfbf-2023a-R-4.3.2
module load SAMtools/1.18-GCC-12.3.0
module load BCFtools/1.18-GCC-12.3.0

WORKDIR=/scratch/eeg37520/doliolid_popgen/mitobim/vcf_dgeg
REF=${WORKDIR}/dgeg_ref.fasta

cd ${WORKDIR}

# ---- Verify all BAMs exist ----
echo "=== Checking BAMs ==="
MISSING=0
while read SAMP; do
    BAM=${WORKDIR}/bams/${SAMP}.dgeg.markdup.bam
    if [[ ! -f ${BAM} ]]; then
        echo "MISSING: ${BAM}"
        MISSING=$((MISSING + 1))
    fi
done < ${WORKDIR}/sample_names_105.txt

if [[ ${MISSING} -gt 0 ]]; then
    echo "ERROR: ${MISSING} BAMs missing. Check remap array job logs."
    exit 1
fi

# ---- Create BAM list ----
while read SAMP; do
    echo "${WORKDIR}/bams/${SAMP}.dgeg.markdup.bam"
done < ${WORKDIR}/sample_names_105.txt > ${WORKDIR}/bamlist_105_dgeg.txt

# ---- Joint variant calling ----
echo "=== Running freebayes ==="
freebayes \
    -f ${REF} \
    -L ${WORKDIR}/bamlist_105_dgeg.txt \
    --ploidy 1 \
    --min-mapping-quality 20 \
    --min-base-quality 20 \
    --min-coverage 3 \
    > ${WORKDIR}/mt_105_dgeg_raw.vcf

echo "Raw variants: $(grep -cv '^#' ${WORKDIR}/mt_105_dgeg_raw.vcf)"

# ---- Filter VCF ----
echo "=== Filtering ==="
bcftools filter \
    -e 'QUAL<20' \
    ${WORKDIR}/mt_105_dgeg_raw.vcf \
| bcftools view \
    -m2 -M2 \
    --types snps \
    -o ${WORKDIR}/mt_105_dgeg_filtered.vcf

echo "Filtered SNPs: $(grep -cv '^#' ${WORKDIR}/mt_105_dgeg_filtered.vcf)"

# ---- Extract SNP positions ----
grep -v '^#' ${WORKDIR}/mt_105_dgeg_filtered.vcf \
    | awk '{print $2}' \
    > ${WORKDIR}/snp_positions_dgeg_105.txt

echo "SNP positions written to: ${WORKDIR}/snp_positions_dgeg_105.txt"

# ---- Summary ----
echo ""
echo "=== VCF SUMMARY ==="
echo "Reference: DD_21_05 (D. gegenbauri, 15329 bp)"
echo "Samples: 105"
echo "Raw variants: $(grep -cv '^#' ${WORKDIR}/mt_105_dgeg_raw.vcf)"
echo "Filtered SNPs: $(grep -cv '^#' ${WORKDIR}/mt_105_dgeg_filtered.vcf)"
echo "Output: ${WORKDIR}/mt_105_dgeg_filtered.vcf"
echo "=== Done ==="
