"""Tests for mem_matesw functionality and paired-end alignment."""

import pytest
from pathlib import Path
from bwamem import BwaAligner

# Test sequences (already MK converted + RC for read2)
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


@pytest.fixture
def pe_results(aligner):
    """Run paired-end alignment and return results."""
    # Returns list of (h1, h2, is_proper, isize)
    return aligner.align(SEQ1_CONV, SEQ2_CONV)


def test_matesw_produces_alignments(pe_results):
    """Test that mem_matesw produces alignment results."""
    assert len(pe_results) > 0, "Should produce at least one paired alignment"


def test_matesw_rRNA_hits(pe_results):
    """Test that mate rescue finds rRNA hits."""
    # hit is (ctg, r_st, r_en, strand, q_st, q_en, mapq, cigar_str, NM, score)
    read1_rRNA = sum(1 for pa in pe_results if pa[0] and "rRNA" in pa[0][0])
    read2_rRNA = sum(1 for pa in pe_results if pa[1] and "rRNA" in pa[1][0])

    assert read1_rRNA > 0, "Read1 should have rRNA hits"
    assert read2_rRNA > 0, "Read2 should have rRNA hits"


def test_matesw_alignment_structure(pe_results):
    """Test that alignments have the expected structure."""
    # Check first alignment
    first_aln = pe_results[0]
    assert first_aln[0] is not None or first_aln[1] is not None, (
        "Alignment should have at least one read"
    )

    # Check alignment properties
    if first_aln[0]:
        assert isinstance(first_aln[0][0], str), "Read1 should have contig name"
        assert isinstance(first_aln[0][9], int), "Read1 should have score"

    if first_aln[1]:
        assert isinstance(first_aln[1][0], str), "Read2 should have contig name"
        assert isinstance(first_aln[1][9], int), "Read2 should have score"


def test_matesw_same_contig_pairs(pe_results):
    """Test that some pairs map to the same contig."""
    same_contig_pairs = [
        pa for pa in pe_results if pa[0] and pa[1] and pa[0][0] == pa[1][0]
    ]

    assert len(same_contig_pairs) > 0, (
        "Should have some pairs mapping to the same contig"
    )

    # Check that these pairs have valid coordinate relationships
    for pa in same_contig_pairs:
        # pa[0] is h1, pa[1] is h2
        r1_start, r1_end = pa[0][1], pa[0][2]
        r2_start, r2_end = pa[1][1], pa[1][2]
        assert r1_end > r1_start or r2_end > r2_start, (
            "Pairs should have valid coordinate ranges"
        )


def test_matesw_rRNA_top_scores(pe_results):
    """Test that rRNA alignments have reasonable scores."""
    rRNA_pairs = [pa for pa in pe_results if pa[0] and "rRNA" in pa[0][0]]

    assert len(rRNA_pairs) > 0, "Should have rRNA pairs"

    for pa in rRNA_pairs:
        if pa[0]:
            assert pa[0][9] > 0, "Read1 rRNA alignment should have positive score"
        if pa[1]:
            assert pa[1][9] > 0, "Read2 alignment should have positive score"
