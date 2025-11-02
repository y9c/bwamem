#!/usr/bin/env python3
"""
Direct test of mem_matesw in bwamem without coralsnake.

NOTE: This is a standalone script for manual testing/debugging.
For automated testing, see test_matesw.py and test_visualization.py
"""

from pathlib import Path
from bwamem import BwaAligner, visualize_paired_alignment

# Converted sequences (already MK converted + RC for read2)
seq1_conv = "TTTTGGTTTTGGGTGGGGGTTGTTGGGGGGGGTGTTGTGGGGGTGGTT"
seq2_conv = "GGTTTGTGGTGGTTGGGTGTTTGTGGTGGTGTTGTTTTTTGGTTTTTTG"

# Use subset reference in tests folder
test_dir = Path(__file__).parent
mk_prefix = str(test_dir / "test_data" / "reference" / "ref.mk.subset")

print("="*80)
print("DIRECT BWAMEM TEST - NO CORALSNAKE")
print("="*80)
print(f"seq1_conv: {seq1_conv}")
print(f"seq2_conv: {seq2_conv}")
print(f"index: {mk_prefix}")
print("="*80)

aligner = BwaAligner(mk_prefix, min_seed_len=14, max_occ=1000, min_score=20)

print("\nCalling aligner.align(seq1_conv, seq2_conv)...")
try:
    pe_alignments = tuple(aligner.align(seq1_conv, seq2_conv))
    print(f"\nResult: {len(pe_alignments)} paired alignments")
except Exception as e:
    print(f"\n✗ Exception during alignment: {e}")
    import traceback
    traceback.print_exc()
    pe_alignments = []

# Check for rRNA hits
read1_rRNA = sum(1 for pa in pe_alignments if pa.read1 and "rRNA" in pa.read1.ctg)
read2_rRNA = sum(1 for pa in pe_alignments if pa.read2 and "rRNA" in pa.read2.ctg)

print(f"Read1 rRNA hits: {read1_rRNA}")
print(f"Read2 rRNA hits: {read2_rRNA}")

if read1_rRNA == 0:
    print("\n✗ MATE RESCUE FAILED: Read1 has no rRNA hits")
    print("Top 5 Read1 hits:")
    seen = set()
    for pa in pe_alignments[:100]:
        if pa.read1:
            key = (pa.read1.ctg, pa.read1.score)
            if key not in seen:
                seen.add(key)
                print(f"  {pa.read1.ctg}: score={pa.read1.score}")
                if len(seen) >= 5:
                    break
else:
    print("\n✓ SUCCESS: Read1 has rRNA hits!")

# Visualize a few top alignments
print("\n" + "="*80)
print("VISUALIZING TOP ALIGNMENTS")
print("="*80)

# Find rRNA hits for visualization
rRNA_pairs = [pa for pa in pe_alignments if pa.read1 and "rRNA" in pa.read1.ctg]

if rRNA_pairs:
    print(f"\nFound {len(rRNA_pairs)} pairs with rRNA hits. Visualizing top 3:\n")
    for i, pa in enumerate(rRNA_pairs[:3], 1):
        print(f"\n{'='*80}")
        print(f"ALIGNMENT #{i} (Score: Read1={pa.read1.score if pa.read1 else 'N/A'}, "
              f"Read2={pa.read2.score if pa.read2 else 'N/A'})")
        print(f"{'='*80}")
        print(visualize_paired_alignment(pa, seq1_conv, seq2_conv, aligner, line_width=80))
        print()
else:
    # Visualize top 3 alignments if no rRNA hits
    print(f"\nVisualizing top 3 alignments:\n")
    for i, pa in enumerate(pe_alignments[:3], 1):
        print(f"\n{'='*80}")
        print(f"ALIGNMENT #{i} (Score: Read1={pa.read1.score if pa.read1 else 'N/A'}, "
              f"Read2={pa.read2.score if pa.read2 else 'N/A'})")
        print(f"{'='*80}")
        print(visualize_paired_alignment(pa, seq1_conv, seq2_conv, aligner, line_width=80))
        print()

print("\n" + "="*80)
print("Check /tmp/mem_matesw_debug.txt for C debug output")
print("="*80)

