#!/usr/bin/env python3

from bwamem import BwaAligner, BwaIndexer

indexer = BwaIndexer()
# Available algorithms: {['auto', 'rb2', 'bwtsw', 'is']}
# Default block size: {indexer.block_size}

# Example of how to build an index (would need actual FASTA file)
from pathlib import Path

test_dir = Path(__file__).parent.parent
rRNA_fa = test_dir / "test_data" / "rRNA" / "rRNA.fa"
indexer = BwaIndexer(algorithm="auto")
index_path = indexer.build_index(
    str(rRNA_fa), prefix=str(test_dir / "test_data" / "rRNA" / "rRNA")
)

aligner = BwaAligner(index_path, options="-T 0 -k 10", insert_model=(150, 15, 1000))

seq1 = "CCTGTCGCTGGAGAGGTTGGGCCT"
seq2 = "GACGCGCGAGAGAACAGCAGGCCCGC"
# Enforce strict pairing with an explicit insert-size window
pe_alignments = aligner.align(seq1, seq2)


def visualize_pair(aln1, aln2, len1, len2, width=80, pad=10):
    start = min(aln1.pos, aln2.pos)
    end = max(aln1.pos + len1, aln2.pos + len2)
    start = max(0, start - pad)
    end = end + pad
    span = max(1, end - start)

    def place(pos, orient, label):
        idx = int((pos - start) * width / span)
        idx = max(0, min(width - 1, idx))
        line = ["-"] * width
        line[idx] = ">" if orient == "+" else "<"
        s = "".join(line)
        return f"{label}:|{s}|"

    print(f"  span: {start}..{end} (~{span}bp)")
    print("  " + place(aln1.pos, aln1.orient, "R1"))
    print("  " + place(aln2.pos, aln2.orient, "R2"))


for pe_aln in filter(lambda p: p.is_proper_pair, pe_alignments):
    r1 = pe_aln.read1
    r2 = pe_aln.read2
    print(f"  Read1: {r1.rname}:{r1.pos} {r1.orient} {r1.cigar}")
    print(f"  Read2: {r2.rname}:{r2.pos} {r2.orient} {r2.cigar}")
    print(f"  Proper pair: {pe_aln.is_proper_pair}, Insert size: {pe_aln.insert_size}")
    visualize_pair(r1, r2, len(seq1), len(seq2))
