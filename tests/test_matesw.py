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
def paired_alignments(aligner):
    """Run paired-end alignment and return results."""
    return tuple(aligner.align(SEQ1_CONV, SEQ2_CONV))


def test_matesw_produces_alignments(paired_alignments):
    """Test that mem_matesw produces alignment results."""
    assert len(paired_alignments) > 0, "Should produce at least one paired alignment"
    assert len(paired_alignments) >= 700, f"Expected at least 700 alignments, got {len(paired_alignments)}"


def test_matesw_rRNA_hits(paired_alignments):
    """Test that mate rescue finds rRNA hits for read1."""
    read1_rRNA = sum(1 for pa in paired_alignments if pa.read1 and "rRNA" in pa.read1.ctg)
    read2_rRNA = sum(1 for pa in paired_alignments if pa.read2 and "rRNA" in pa.read2.ctg)
    
    assert read1_rRNA > 0, "Read1 should have rRNA hits (mate rescue should work)"
    assert read2_rRNA > 0, "Read2 should have rRNA hits"
    
    # Verify we have a reasonable number of rRNA hits
    assert read1_rRNA >= 10, f"Expected at least 10 rRNA hits for read1, got {read1_rRNA}"


def test_matesw_alignment_structure(paired_alignments):
    """Test that alignments have the expected structure."""
    # Check first alignment
    first_aln = paired_alignments[0]
    assert first_aln.read1 is not None or first_aln.read2 is not None, "Alignment should have at least one read"
    
    # Check alignment properties
    if first_aln.read1:
        assert hasattr(first_aln.read1, 'ctg'), "Read1 should have contig name"
        assert hasattr(first_aln.read1, 'score'), "Read1 should have score"
        assert hasattr(first_aln.read1, 'r_st'), "Read1 should have reference start"
        assert hasattr(first_aln.read1, 'r_en'), "Read1 should have reference end"
    
    if first_aln.read2:
        assert hasattr(first_aln.read2, 'ctg'), "Read2 should have contig name"
        assert hasattr(first_aln.read2, 'score'), "Read2 should have score"


def test_matesw_primary_alignments(paired_alignments):
    """Test that primary alignments are marked correctly."""
    primary_count = sum(1 for pa in paired_alignments 
                       if (pa.read1 and pa.read1.is_primary) or 
                          (pa.read2 and pa.read2.is_primary))
    
    assert primary_count > 0, "Should have at least some primary alignments"


def test_matesw_same_contig_pairs(paired_alignments):
    """Test that some pairs map to the same contig."""
    same_contig_pairs = [pa for pa in paired_alignments 
                        if pa.read1 and pa.read2 and pa.read1.ctg == pa.read2.ctg]
    
    assert len(same_contig_pairs) > 0, "Should have some pairs mapping to the same contig"
    
    # Check that these pairs have valid coordinate relationships
    for pa in same_contig_pairs[:5]:  # Check first 5
        r1_start = min(pa.read1.r_st, pa.read1.r_en)
        r1_end = max(pa.read1.r_st, pa.read1.r_en)
        r2_start = min(pa.read2.r_st, pa.read2.r_en)
        r2_end = max(pa.read2.r_st, pa.read2.r_en)
        
        # At least one read should have valid coordinates
        assert r1_end > r1_start or r2_end > r2_start, "Pairs should have valid coordinate ranges"


def test_matesw_rRNA_top_scores(paired_alignments):
    """Test that rRNA alignments have reasonable scores."""
    rRNA_pairs = [pa for pa in paired_alignments if pa.read1 and "rRNA" in pa.read1.ctg]
    
    assert len(rRNA_pairs) > 0, "Should have rRNA pairs"
    
    # Check scores are positive and reasonable
    for pa in rRNA_pairs[:10]:  # Check first 10 rRNA pairs
        if pa.read1:
            assert pa.read1.score > 0, "Read1 rRNA alignment should have positive score"
        if pa.read2:
            assert pa.read2.score > 0, "Read2 alignment should have positive score"

