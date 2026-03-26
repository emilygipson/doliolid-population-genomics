#!/bin/bash
#SBATCH --job-name=hetplat_screen
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Heteroplasmy screen: remap 5 high-coverage samples to D. gegenbauri
# reference (DD_21_05) and call per-site allele frequencies with bcftools
#
# Samples: DD_11_24, DD_11_11, DD_9_11, DD_11_05, DD_9_06
# Detection threshold: minor allele freq >= 5%, min depth 20x

module load BWA/0.7.18-GCCcore-12.3.0
module load SAMtools/1.18-GCC-12.3.0
module load BCFtools/1.18-GCC-12.3.0

WORKDIR=/scratch/eeg37520/doliolid_popgen/mitobim/hetplat_screen
READS=/scratch/eeg37520/doliolid_popgen/mitobim/mt_reads/fastq
REF=${WORKDIR}/dgeg_ref.fasta

SAMPLES="DD_11_24 DD_11_11 DD_9_11 DD_11_05 DD_9_06"

cd ${WORKDIR}

# ---- Index reference ----
bwa index ${REF}
samtools faidx ${REF}

# ---- Remap each sample ----
for SAMP in ${SAMPLES}; do
    echo "=== Processing ${SAMP} ==="

    R1=${READS}/${SAMP}_R1.fastq.gz
    R2=${READS}/${SAMP}_R2.fastq.gz

    if [[ ! -f ${R1} || ! -f ${R2} ]]; then
        echo "WARNING: reads not found for ${SAMP}, skipping"
        continue
    fi

    # Map to D. gegenbauri reference
    bwa mem -t 8 \
        -R "@RG\tID:${SAMP}\tSM:${SAMP}\tPL:ILLUMINA" \
        ${REF} ${R1} ${R2} \
    | samtools fixmate -m - - \
    | samtools sort -@ 4 -o ${WORKDIR}/${SAMP}.dgeg.sorted.bam

    # Mark duplicates
    samtools markdup ${WORKDIR}/${SAMP}.dgeg.sorted.bam \
        ${WORKDIR}/${SAMP}.dgeg.markdup.bam
    samtools index ${WORKDIR}/${SAMP}.dgeg.markdup.bam

    # Coverage stats
    samtools coverage ${WORKDIR}/${SAMP}.dgeg.markdup.bam \
        > ${WORKDIR}/${SAMP}.dgeg.coverage.txt

    # Clean intermediate
    rm ${WORKDIR}/${SAMP}.dgeg.sorted.bam

    echo "=== ${SAMP} done ==="
done

# ---- Call allele frequencies per sample ----
echo "=== Calling allele frequencies ==="

for SAMP in ${SAMPLES}; do
    BAM=${WORKDIR}/${SAMP}.dgeg.markdup.bam
    if [[ ! -f ${BAM} ]]; then continue; fi

    # mpileup with per-base quality, then call with ploidy 1 to get all variants
    # Use -d 0 for no depth cap, -q 20 min mapping qual, -Q 20 min base qual
    bcftools mpileup \
        -f ${REF} \
        -d 0 \
        -q 20 \
        -Q 20 \
        --annotate FORMAT/AD,FORMAT/DP \
        ${BAM} \
    | bcftools call \
        --ploidy 1 \
        -m \
        -Ov \
        -o ${WORKDIR}/${SAMP}.dgeg.raw.vcf

    # Also generate per-site depth and allele counts for heteroplasmy parsing
    bcftools mpileup \
        -f ${REF} \
        -d 0 \
        -q 20 \
        -Q 20 \
        --annotate FORMAT/AD,FORMAT/DP \
        ${BAM} \
    | bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP\t%AD]\n' \
        > ${WORKDIR}/${SAMP}.dgeg.allele_depths.tsv

done

# ---- Screen for heteroplasmic sites ----
echo "=== Screening for heteroplasmy ==="
echo -e "Sample\tPos\tRef\tAlt\tDepth\tRef_count\tAlt_count\tMinor_AF" \
    > ${WORKDIR}/heteroplasmy_candidates.tsv

for SAMP in ${SAMPLES}; do
    ADFILE=${WORKDIR}/${SAMP}.dgeg.allele_depths.tsv
    if [[ ! -f ${ADFILE} ]]; then continue; fi

    awk -v samp="${SAMP}" 'BEGIN{OFS="\t"} {
        # Parse AD field (ref,alt counts)
        split($6, ad, ",")
        ref_ct = ad[1]
        alt_ct = 0
        for (i=2; i<=length(ad); i++) alt_ct += ad[i]
        total = ref_ct + alt_ct
        if (total < 20) next
        if (alt_ct == 0) next
        minor = (alt_ct < ref_ct) ? alt_ct : ref_ct
        maf = minor / total
        if (maf >= 0.05)
            print samp, $2, $3, $4, total, ref_ct, alt_ct, sprintf("%.3f", maf)
    }' ${ADFILE} >> ${WORKDIR}/heteroplasmy_candidates.tsv

done

# ---- Summary ----
echo ""
echo "=== HETEROPLASMY SCREEN SUMMARY ==="
echo "Samples screened: $(echo ${SAMPLES} | wc -w)"
echo "Detection threshold: MAF >= 5%, depth >= 20x"
echo ""

# Count per sample
for SAMP in ${SAMPLES}; do
    n=$(grep -c "^${SAMP}" ${WORKDIR}/heteroplasmy_candidates.tsv 2>/dev/null || echo 0)
    cov=$(awk '{print $4}' ${WORKDIR}/${SAMP}.dgeg.coverage.txt 2>/dev/null | tail -1)
    echo "${SAMP}: ${n} candidate heteroplasmic sites (mean cov: ${cov}x)"
done

echo ""
echo "Total candidates: $(tail -n+2 ${WORKDIR}/heteroplasmy_candidates.tsv | wc -l)"
echo "Results in: ${WORKDIR}/heteroplasmy_candidates.tsv"
echo "=== Done ==="
