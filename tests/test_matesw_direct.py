"""Direct test of mem_matesw with detailed output and visualization.

This test provides detailed output similar to the original manual test script,
but runs as part of the pytest suite. Useful for debugging and verification.
"""

import pytest
from pathlib import Path
from bwamem import BwaAligner, visualize_paired_alignment

# Converted sequences (already MK converted + RC for read2)
SEQ1_CONV = "TTTTGGTTTTGGGTGGGGGTTGTTGGGGGGGGTGTTGTGGGGGTGGTT"
SEQ2_CONV = "GGTTTGTGGTGGTTGGGTGTTTGTGGTGGTGTTGTTTTTTGGTTTTTTG"


@pytest.fixture
def test_index_path():
    """Return path to test reference index."""
    test_dir = Path(__file__).parent
    return str(test_dir / "test_data" / "reference" / "ref.mk.subset")


@pytest.fixture
def aligner(test_index_path):
    """Create BwaAligner instance for testing."""
    return BwaAligner(test_index_path, min_seed_len=14, max_occ=1000, min_score=20)


def test_matesw_direct_with_output(aligner):
    """Test mem_matesw with detailed output - shows alignment process."""
    print("="*80)
    print("DIRECT BWAMEM TEST - NO CORALSNAKE")
    print("="*80)
    print(f"seq1_conv: {SEQ1_CONV}")
    print(f"seq2_conv: {SEQ2_CONV}")
    print(f"index: {aligner}")
    print("="*80)

    print("\nCalling aligner.align(seq1_conv, seq2_conv)...")
    pe_alignments = tuple(aligner.align(SEQ1_CONV, SEQ2_CONV))
    print(f"\nResult: {len(pe_alignments)} paired alignments")
    
    assert len(pe_alignments) > 0, "Should produce at least one paired alignment"


def test_matesw_direct_rRNA_hits(aligner):
    """Test mem_matesw rRNA hits with detailed output."""
    pe_alignments = tuple(aligner.align(SEQ1_CONV, SEQ2_CONV))
    
    # Check for rRNA hits
    read1_rRNA = sum(1 for pa in pe_alignments if pa.read1 and "rRNA" in pa.read1.ctg)
    read2_rRNA = sum(1 for pa in pe_alignments if pa.read2 and "rRNA" in pa.read2.ctg)

    print(f"\nRead1 rRNA hits: {read1_rRNA}")
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
    
    # Assertions for pytest
    assert read1_rRNA > 0, "Read1 should have rRNA hits (mate rescue should work)"
    assert read2_rRNA > 0, "Read2 should have rRNA hits"


def test_matesw_direct_visualization(aligner):
    """Test visualization of top alignments with detailed output."""
    pe_alignments = tuple(aligner.align(SEQ1_CONV, SEQ2_CONV))
    
    print("\n" + "="*80)
    print("VISUALIZING TOP ALIGNMENTS")
    print("="*80)

    # Find rRNA hits for visualization
    rRNA_pairs = [pa for pa in pe_alignments if pa.read1 and "rRNA" in pa.read1.ctg]

    if rRNA_pairs:
        print(f"\nFound {len(rRNA_pairs)} pairs with rRNA hits. Visualizing top 3:\n")
        visualize_pairs = rRNA_pairs[:3]
    else:
        # Visualize top 3 alignments if no rRNA hits
        print(f"\nVisualizing top 3 alignments:\n")
        visualize_pairs = pe_alignments[:3]
    
    assert len(visualize_pairs) > 0, "Should have at least one alignment to visualize"
    
    for i, pa in enumerate(visualize_pairs, 1):
        print(f"\n{'='*80}")
        print(f"ALIGNMENT #{i} (Score: Read1={pa.read1.score if pa.read1 else 'N/A'}, "
              f"Read2={pa.read2.score if pa.read2 else 'N/A'})")
        print(f"{'='*80}")
        try:
            viz_output = visualize_paired_alignment(pa, SEQ1_CONV, SEQ2_CONV, aligner, line_width=80)
            print(viz_output)
            assert len(viz_output) > 0, "Visualization should produce output"
        except Exception as e:
            print(f"Visualization failed: {e}")
            import traceback
            traceback.print_exc()
            pytest.fail(f"Visualization failed: {e}")
        print()

    print("\n" + "="*80)
    print("Check /tmp/mem_matesw_debug.txt for C debug output")
    print("="*80)

