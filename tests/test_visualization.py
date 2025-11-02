"""Tests for alignment visualization functionality."""

import pytest
from pathlib import Path
from bwamem import BwaAligner, visualize_paired_alignment, Alignment

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


@pytest.fixture
def rRNA_pair(paired_alignments):
    """Get a paired alignment with rRNA hits."""
    rRNA_pairs = [pa for pa in paired_alignments if pa.read1 and "rRNA" in pa.read1.ctg]
    assert len(rRNA_pairs) > 0, "Should have at least one rRNA pair"
    return rRNA_pairs[0]


def test_alignment_visualize_method(rRNA_pair, aligner):
    """Test that Alignment.visualize() method works correctly."""
    read1 = rRNA_pair.read1
    assert read1 is not None, "Read1 should exist"
    
    # Helper function to get reference sequence
    def ref_seq_fn(ctg, start, end):
        return aligner.seq(ctg, start, end)
    
    # Test visualization
    viz_output = read1.visualize(SEQ1_CONV, ref_seq_fn, line_width=80)
    
    # Check output structure
    assert isinstance(viz_output, str), "Visualization should return a string"
    assert len(viz_output) > 0, "Visualization should not be empty"
    
    # Check that it contains expected elements
    assert read1.ctg in viz_output, "Visualization should contain contig name"
    assert "Ref" in viz_output, "Visualization should contain reference line"
    assert "Qry" in viz_output, "Visualization should contain query line"
    
    # Check for match indicators
    assert "|" in viz_output or "X" in viz_output, "Visualization should contain match indicators"


def test_alignment_visualize_formatting(rRNA_pair, aligner):
    """Test that visualization has correct formatting and alignment."""
    read1 = rRNA_pair.read1
    
    def ref_seq_fn(ctg, start, end):
        return aligner.seq(ctg, start, end)
    
    viz_output = read1.visualize(SEQ1_CONV, ref_seq_fn, line_width=80)
    lines = viz_output.split('\n')
    
    # Find reference and match lines
    ref_line = None
    match_line = None
    query_line = None
    
    for line in lines:
        if line.startswith('Ref  '):
            ref_line = line
        elif line.strip().startswith('|') or line.strip().startswith('X'):
            match_line = line
        elif line.startswith('Qry  '):
            query_line = line
    
    assert ref_line is not None, "Should have reference line"
    assert match_line is not None, "Should have match line"
    assert query_line is not None, "Should have query line"
    
    # Check that match line is properly aligned with reference sequence
    # The match indicators should start at the same position as the sequence in ref_line
    ref_seq_start = ref_line.find('G', ref_line.find(' '))  # Find first G after space
    match_indicator_start = len(match_line) - len(match_line.lstrip())
    
    assert ref_seq_start == match_indicator_start, \
        f"Match indicators should align with sequence start (ref starts at {ref_seq_start}, match at {match_indicator_start})"


def test_visualize_paired_alignment_function(rRNA_pair, aligner):
    """Test visualize_paired_alignment() function."""
    viz_output = visualize_paired_alignment(rRNA_pair, SEQ1_CONV, SEQ2_CONV, aligner, line_width=80)
    
    # Check output structure
    assert isinstance(viz_output, str), "Visualization should return a string"
    assert len(viz_output) > 0, "Visualization should not be empty"
    
    # Check that it contains expected sections
    assert "PAIRED-END ALIGNMENT VISUALIZATION" in viz_output
    assert "--- READ 1 ---" in viz_output or "--- READ 2 ---" in viz_output
    
    # Check for both reads if they exist
    if rRNA_pair.read1:
        assert "--- READ 1 ---" in viz_output, "Should visualize Read1"
    if rRNA_pair.read2:
        assert "--- READ 2 ---" in viz_output, "Should visualize Read2"


def test_visualize_paired_alignment_content(rRNA_pair, aligner):
    """Test that paired alignment visualization contains expected content."""
    viz_output = visualize_paired_alignment(rRNA_pair, SEQ1_CONV, SEQ2_CONV, aligner)
    
    lines = viz_output.split('\n')
    
    # Check for proper pair information
    has_proper_pair_info = any("Proper pair" in line for line in lines)
    assert has_proper_pair_info, "Should contain proper pair information"
    
    # Check that both reads are visualized if they exist
    read1_section = False
    read2_section = False
    
    in_read1 = False
    in_read2 = False
    for line in lines:
        if "--- READ 1 ---" in line:
            in_read1 = True
            read1_section = True
        if "--- READ 2 ---" in line:
            in_read2 = True
            read2_section = True
            in_read1 = False
        
        if in_read1 and line.startswith('Ref  '):
            assert "Ref" in line and len(line) > 10, "Read1 should have reference line"
        if in_read2 and line.startswith('Ref  '):
            assert "Ref" in line and len(line) > 10, "Read2 should have reference line"
    
    if rRNA_pair.read1:
        assert read1_section, "Should visualize Read1 section"
    if rRNA_pair.read2:
        assert read2_section, "Should visualize Read2 section"


def test_visualize_paired_alignment_pairing_info(rRNA_pair, aligner):
    """Test that pairing information is displayed correctly."""
    viz_output = visualize_paired_alignment(rRNA_pair, SEQ1_CONV, SEQ2_CONV, aligner)
    
    # If both reads are on the same contig, should show pairing info
    if rRNA_pair.read1 and rRNA_pair.read2:
        if rRNA_pair.read1.ctg == rRNA_pair.read2.ctg:
            assert "--- PAIRING INFO ---" in viz_output, "Should show pairing info for same-contig pairs"
            assert "Gap between reads" in viz_output, "Should show gap information"
            assert "Insert size" in viz_output, "Should show insert size information"


def test_alignment_reverse_complement_handling(paired_alignments, aligner):
    """Test that reverse complement sequences are handled correctly in visualization."""
    # Find an alignment with reverse strand
    reverse_aln = None
    for pa in paired_alignments:
        if pa.read1 and pa.read1.strand < 0:
            reverse_aln = pa.read1
            break
    
    if reverse_aln:
        def ref_seq_fn(ctg, start, end):
            return aligner.seq(ctg, start, end)
        
        viz_output = reverse_aln.visualize(SEQ1_CONV, ref_seq_fn)
        # Should handle reverse complement correctly
        assert len(viz_output) > 0, "Reverse complement alignment should visualize"


def test_visualization_multiple_alignments(rRNA_pair, aligner):
    """Test visualizing multiple alignments."""
    # Visualize first rRNA pair
    viz1 = visualize_paired_alignment(rRNA_pair, SEQ1_CONV, SEQ2_CONV, aligner)
    assert len(viz1) > 0
    
    # Should be able to visualize multiple times
    viz2 = visualize_paired_alignment(rRNA_pair, SEQ1_CONV, SEQ2_CONV, aligner)
    assert len(viz2) > 0
    assert len(viz1) == len(viz2), "Multiple visualizations should be consistent"


def test_visualization_line_width_parameter(rRNA_pair, aligner):
    """Test that line_width parameter works correctly."""
    # Test with different line widths
    for line_width in [40, 60, 80, 100]:
        viz_output = visualize_paired_alignment(rRNA_pair, SEQ1_CONV, SEQ2_CONV, aligner, line_width=line_width)
        assert len(viz_output) > 0, f"Visualization should work with line_width={line_width}"
        
        # Check that lines don't exceed line_width (allowing for labels)
        lines = viz_output.split('\n')
        for line in lines:
            if line.startswith('Ref  ') or line.startswith('Qry  '):
                # Allow for position labels and some overflow
                assert len(line) <= line_width + 30, f"Line should respect line_width approximately: {len(line)}"

