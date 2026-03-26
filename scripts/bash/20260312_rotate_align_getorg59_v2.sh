#!/bin/bash
#SBATCH --job-name=rotate_getorg59_v2
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=/scratch/eeg37520/doliolid_popgen/mitobim/getorg_complete/logs/rotate_align_v2_%j.out
#SBATCH --error=/scratch/eeg37520/doliolid_popgen/mitobim/getorg_complete/logs/rotate_align_v2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu

module load MAFFT/7.520-GCC-12.3.0-with-extensions

GETORG=/scratch/eeg37520/doliolid_popgen/mitobim/getorg_complete
OUTDIR=$GETORG/rotated_alignment
mkdir -p $OUTDIR

echo "[$(date)] Rotating assemblies to coding region anchor..."

python3 << 'PYEOF'
import os
import sys
import glob

# Coding region anchor from DD_21_05 position 5000-5030
# (conserved across samples, not in AT-rich control region)
anchor = "TAGTGACTCAGCCTACAGGTCCGTTTTCTT"
getorg_dir = "/scratch/eeg37520/doliolid_popgen/mitobim/getorg_complete"
out_dir = os.path.join(getorg_dir, "rotated_alignment")

def read_fasta(path):
    seq = ""
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq.upper()

def rc(seq):
    comp = str.maketrans("ACGTNRYSWKMBDHV", "TGCANYRSWMKVHDB")
    return seq.translate(comp)[::-1]

def find_anchor(seq, anchor, max_mm=5):
    """Search for anchor allowing up to max_mm mismatches, try fwd then rc"""
    # Exact forward
    pos = seq.find(anchor)
    if pos >= 0:
        return seq[pos:] + seq[:pos], "fwd", pos, 0
    
    # Exact reverse complement
    rc_seq = rc(seq)
    pos = rc_seq.find(anchor)
    if pos >= 0:
        return rc_seq[pos:] + rc_seq[:pos], "rc", pos, 0
    
    # Fuzzy forward
    best_fwd = (None, len(anchor))  # (pos, mismatches)
    for i in range(len(seq) - len(anchor)):
        mm = sum(1 for a, b in zip(seq[i:i+len(anchor)], anchor) if a != b)
        if mm < best_fwd[1]:
            best_fwd = (i, mm)
    
    # Fuzzy rc
    best_rc = (None, len(anchor))
    for i in range(len(rc_seq) - len(anchor)):
        mm = sum(1 for a, b in zip(rc_seq[i:i+len(anchor)], anchor) if a != b)
        if mm < best_rc[1]:
            best_rc = (i, mm)
    
    # Pick best match
    if best_fwd[1] <= best_rc[1] and best_fwd[1] <= max_mm:
        pos = best_fwd[0]
        return seq[pos:] + seq[:pos], f"fwd_fuzzy({best_fwd[1]}mm)", pos, best_fwd[1]
    elif best_rc[1] <= max_mm:
        pos = best_rc[0]
        return rc_seq[pos:] + rc_seq[:pos], f"rc_fuzzy({best_rc[1]}mm)", pos, best_rc[1]
    
    return seq, "no_match", -1, min(best_fwd[1], best_rc[1])

fastas = sorted(glob.glob(os.path.join(getorg_dir, "*_mt.fasta")))
print(f"Found {len(fastas)} assemblies", file=sys.stderr)

rotated_file = os.path.join(out_dir, "all_59_rotated_unaligned.fasta")
n_success = 0
n_fail = 0

with open(rotated_file, "w") as out:
    for fa in fastas:
        sample = os.path.basename(fa).replace("_mt.fasta", "")
        seq = read_fasta(fa)
        rotated, orient, pos, mm = find_anchor(seq, anchor, max_mm=5)
        print(f"  {sample}: {orient} at {pos}, len={len(seq)}, mm={mm}", file=sys.stderr)
        if orient == "no_match":
            print(f"  WARNING: {sample} not rotated (best match had {mm} mismatches)", file=sys.stderr)
            n_fail += 1
        else:
            n_success += 1
        out.write(f">{sample}\n{rotated}\n")

print(f"\nRotated: {n_success}, Failed: {n_fail}", file=sys.stderr)
print(f"Wrote to {rotated_file}", file=sys.stderr)
PYEOF

# Align with MAFFT
echo ""
echo "[$(date)] Aligning with MAFFT..."
mafft --auto --thread 8 $OUTDIR/all_59_rotated_unaligned.fasta > $OUTDIR/all_59_getorg_rotated_aligned.fasta

# Alignment stats
echo ""
echo "[$(date)] Alignment stats:"
NSEQ=$(grep -c ">" $OUTDIR/all_59_getorg_rotated_aligned.fasta)
ALEN=$(awk '/^>/{if(l) print l; l=0; next} {l+=length($0)} END{print l}' $OUTDIR/all_59_getorg_rotated_aligned.fasta | sort -u)
echo "  Sequences: $NSEQ"
echo "  Alignment length: $ALEN"
echo ""
echo "[$(date)] Done."
