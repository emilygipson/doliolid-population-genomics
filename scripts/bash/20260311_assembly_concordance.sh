#!/bin/bash
#SBATCH --job-name=asm_concordance
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/logs/concordance_%j.out
#SBATCH --error=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/logs/concordance_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu

module load MAFFT/7.520-GCC-12.3.0-with-extensions

OVERLAP=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/overlap_samples.txt
GETORG=/scratch/eeg37520/doliolid_popgen/mitobim/getorg_complete
MITOBIM=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/assemblies
OUTDIR=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/concordance
mkdir -p $OUTDIR

echo "sample	getorg_len	mitobim_len	aligned_len	matches	mismatches	gaps	pct_identity" > $OUTDIR/concordance_results.tsv

while read SAMPLE; do
  echo "[$(date)] Processing $SAMPLE..."

  # Get GetOrganelle assembly (single contig, strip header)
  GETORG_FA=$GETORG/${SAMPLE}_mt.fasta
  
  # Get MITObim conspecific assembly (noIUPAC file)
  MITOBIM_FA=$(find $MITOBIM/$SAMPLE -name "*noIUPAC.fasta" | head -1)

  if [[ ! -f "$GETORG_FA" || ! -f "$MITOBIM_FA" ]]; then
    echo "  MISSING: $SAMPLE"
    continue
  fi

  # Create temp file with both sequences, relabeled
  TMPFILE=$OUTDIR/${SAMPLE}_pair.fasta
  echo ">${SAMPLE}_getorg" > $TMPFILE
  grep -v "^>" $GETORG_FA | tr -d '\n' >> $TMPFILE
  echo "" >> $TMPFILE
  echo ">${SAMPLE}_mitobim" >> $TMPFILE
  grep -v "^>" $MITOBIM_FA | tr -d '\n' >> $TMPFILE
  echo "" >> $TMPFILE

  # Get raw lengths
  GLEN=$(grep -v "^>" $GETORG_FA | tr -d '\n' | wc -c)
  MLEN=$(grep -v "^>" $MITOBIM_FA | tr -d '\n' | wc -c)

  # Align
  ALNFILE=$OUTDIR/${SAMPLE}_aligned.fasta
  mafft --auto --quiet $TMPFILE > $ALNFILE

  # Count matches, mismatches, gaps using awk on the alignment
  python3 -c "
import sys
seqs = {}
name = None
with open('$ALNFILE') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            name = line[1:]
            seqs[name] = ''
        else:
            seqs[name] += line.upper()
names = list(seqs.keys())
s1, s2 = seqs[names[0]], seqs[names[1]]
matches = mismatches = gaps = 0
for a, b in zip(s1, s2):
    if a == '-' or b == '-':
        gaps += 1
    elif a == b:
        matches += 1
    else:
        mismatches += 1
alen = len(s1)
pct = matches / (matches + mismatches) * 100 if (matches + mismatches) > 0 else 0
print(f'$SAMPLE\t{$GLEN}\t{$MLEN}\t{alen}\t{matches}\t{mismatches}\t{gaps}\t{pct:.4f}')
" >> $OUTDIR/concordance_results.tsv

  # Clean up pair file
  rm $TMPFILE

done < $OVERLAP

echo ""
echo "=========================================="
echo "SUMMARY"
echo "=========================================="
cat $OUTDIR/concordance_results.tsv
echo ""
echo "[$(date)] Done"
