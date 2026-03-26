#!/bin/bash
#SBATCH --job-name=numt_screen
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/logs/numt_screen_%j.out
#SBATCH --error=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/logs/numt_screen_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu

REFSEQ=/scratch/eeg37520/doliolid_popgen/mitobim/getorg_complete/DD_21_05_mt.fasta
ASSEMBLYDIR=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/assemblies
OUTDIR=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/numt_screen
mkdir -p $OUTDIR

# Extract each gene from the reference using coordinates from GFF
python3 << 'PYEOF'
import sys

ref_path = "/scratch/eeg37520/doliolid_popgen/mitobim/getorg_complete/DD_21_05_mt.fasta"
with open(ref_path) as f:
    lines = f.readlines()
ref_seq = "".join(l.strip() for l in lines if not l.startswith(">"))

genes = [
    ("nad4", 229, 1558), ("nad6", 3105, 3591), ("cox3", 4450, 5245),
    ("cox2", 5368, 6034), ("cob", 6017, 7097), ("nad4l", 7341, 7602),
    ("nad3", 7586, 7934), ("nad5", 8454, 10140), ("atp6", 10217, 10847),
    ("cox1", 10929, 12480), ("atp8", 12592, 12739), ("nad1", 12897, 13806),
    ("nad2", 13930, 14938)
]

outdir = "/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/numt_screen"
for gene, start, end in genes:
    gene_seq = ref_seq[start:end]
    with open(f"{outdir}/ref_{gene}.fasta", "w") as f:
        f.write(f">{gene}\n{gene_seq}\n")
    print(f"{gene}: {end-start} bp extracted")

PYEOF

echo "[$(date)] Reference genes extracted"

# Ascidian mitochondrial code (13) translation
# For each sample, extract each gene region and check for premature stops
python3 << 'PYEOF'
import os, glob, sys

codon_table = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"M","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"W","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"G","AGG":"G","GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

genes = [
    ("nad4", 229, 1558), ("nad6", 3105, 3591), ("cox3", 4450, 5245),
    ("cox2", 5368, 6034), ("cob", 6017, 7097), ("nad4l", 7341, 7602),
    ("nad3", 7586, 7934), ("nad5", 8454, 10140), ("atp6", 10217, 10847),
    ("cox1", 10929, 12480), ("atp8", 12592, 12739), ("nad1", 12897, 13806),
    ("nad2", 13930, 14938)
]

# Read reference to get anchor for rotation
ref_path = "/scratch/eeg37520/doliolid_popgen/mitobim/getorg_complete/DD_21_05_mt.fasta"
with open(ref_path) as f:
    lines = f.readlines()
ref_seq = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
anchor = ref_seq[:30]

assembly_dir = "/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/assemblies"
outdir = "/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/numt_screen"

sample_list = "/scratch/eeg37520/doliolid_popgen/mitobim/sample_names_105.txt"
with open(sample_list) as f:
    samples = [l.strip() for l in f if l.strip()]

def read_fasta(path):
    seq = ""
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq.upper()

def rotate_to_anchor(seq, anchor):
    pos = seq.find(anchor)
    if pos >= 0:
        return seq[pos:] + seq[:pos]
    # Try fuzzy match (2 mismatches)
    for i in range(len(seq) - len(anchor)):
        mm = sum(1 for a, b in zip(seq[i:i+len(anchor)], anchor) if a != b)
        if mm <= 2:
            return seq[i:] + seq[:i]
    return None

def translate_check(seq, gene_name):
    """Translate and return list of premature stop codon positions"""
    seq = seq.replace("N", "")
    stops = []
    n_codons = len(seq) // 3
    if n_codons < 2:
        return stops
    for i in range(n_codons - 1):  # skip last codon
        codon = seq[i*3:(i+1)*3]
        if len(codon) == 3 and "N" not in codon:
            aa = codon_table.get(codon, "?")
            if aa == "*":
                stops.append((i+1, codon))
    return stops

# Results
print("sample\tgene\tgene_len\tpremature_stops\tstop_details")

gene_totals = {g[0]: {"clean": 0, "with_stops": 0} for g in genes}

for sample in samples:
    fa_files = glob.glob(f"{assembly_dir}/{sample}/**/noIUPAC.fasta", recursive=True)
    if not fa_files:
        fa_files = glob.glob(f"{assembly_dir}/{sample}/**/*noIUPAC.fasta", recursive=True)
    if not fa_files:
        print(f"{sample}\tALL\t-\t-\tMISSING_ASSEMBLY", file=sys.stderr)
        continue

    raw_seq = read_fasta(fa_files[0])
    rotated = rotate_to_anchor(raw_seq, anchor)
    if rotated is None:
        print(f"{sample}\tALL\t-\t-\tROTATION_FAILED", file=sys.stderr)
        continue

    for gene_name, start, end in genes:
        gene_len = end - start
        gene_seq = rotated[start:end]
        stops = translate_check(gene_seq, gene_name)
        if stops:
            gene_totals[gene_name]["with_stops"] += 1
            details = "; ".join(f"pos{p}:{c}" for p, c in stops)
            print(f"{sample}\t{gene_name}\t{gene_len}\t{len(stops)}\t{details}")
        else:
            gene_totals[gene_name]["clean"] += 1

print("\n========================================")
print("GENE SUMMARY")
print("========================================")
print(f"{'Gene':<8} {'Length':<8} {'Clean':<8} {'With stops':<12}")
total_problems = 0
for gene_name, start, end in genes:
    c = gene_totals[gene_name]["clean"]
    s = gene_totals[gene_name]["with_stops"]
    total_problems += s
    print(f"{gene_name:<8} {end-start:<8} {c:<8} {s:<12}")
print(f"\nTotal samples with any premature stop: {total_problems}")
print(f"Total clean sample-gene combinations: {sum(v['clean'] for v in gene_totals.values())}")

PYEOF

echo "[$(date)] NUMT screen complete"
