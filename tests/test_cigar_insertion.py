"""Test for correct CIGAR generation with insertions.

This test validates that BWA correctly generates CIGAR strings for alignments
with insertions, ensuring that the reference span is calculated correctly.

Bug: Previously, the code would force the reference span to match the query span,
which prevented BWA from correctly detecting insertions via dynamic programming.
This resulted in incorrect CIGARs like "44M" instead of "4M1I39M".
"""

import pytest
import tempfile
import os
from bwamem import BwaIndexer, BwaAligner


def test_cigar_with_insertion():
    """Test that CIGAR correctly represents a 1-base insertion in the query."""
    # Create a synthetic reference with padding before and after the alignment region
    # This ensures BWA properly identifies the insertion without boundary effects
    prefix = 'GGCGAGCCACCGCCCGTCCCCGCCCCTTGCCTCTCGGCGCCCCCTCGATGCTCTTAGCTGAGTGTCCCGCGGGGCCCGAAGCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCCCGAGCCGCCTGGATACCGCAGCTAGGAATAATGGAAT'
    # Alignment region: reference has GAATAGGACC... (43bp without the insertion)
    alignment_region = 'GAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG'
    suffix = 'CCATGATTAAGAGGGACGGCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACCAGAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTC'
    reference = prefix + alignment_region + suffix
    
    # Query has an extra 'A' at position 4 compared to the reference
    # Visual alignment:
    #   Reference: GAAT-AGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG
    #   Query:     GAATAAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG
    #              ^^^^I^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    query = 'GAATAAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG'
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write reference FASTA
        ref_file = os.path.join(tmpdir, 'ref.fa')
        with open(ref_file, 'w') as f:
            f.write('>ref\n')
            f.write(reference + '\n')
        
        # Build index
        indexer = BwaIndexer(verbose=0)
        index_prefix = indexer.build_index(ref_file)
        
        # Align query
        aligner = BwaAligner(index_prefix)
        alignments = aligner.align(query)
        
        # Should find one alignment
        assert len(alignments) == 1, f"Expected 1 alignment, got {len(alignments)}"
        
        aln = alignments[0]
        
        # Check basic alignment properties
        assert aln.ctg == 'ref'
        assert aln.q_st == 0
        assert aln.q_en == 44  # Full query aligned
        # The alignment starts at the position where prefix ends
        expected_start = len(prefix)
        assert aln.r_st == expected_start
        
        # The key assertion: reference span should be 43, not 44
        # because there's a 1-base insertion in the query
        ref_span = aln.r_en - aln.r_st
        assert ref_span == 43, f"Expected reference span of 43, got {ref_span}"
        
        # Check CIGAR string contains insertion operation
        assert 'I' in aln.cigar_str, f"CIGAR should contain insertion: {aln.cigar_str}"
        
        # Check specific CIGAR structure: should be 4M1I39M
        assert aln.cigar_str == '4M1I39M', f"Expected CIGAR '4M1I39M', got '{aln.cigar_str}'"
        
        # Check edit distance is 1 (one insertion)
        assert aln.NM == 1, f"Expected edit distance 1, got {aln.NM}"
        
        # Validate CIGAR consistency
        q_consumed = sum(op_len for op_len, op in aln.cigar if op in [0, 1, 4, 7, 8])
        r_consumed = sum(op_len for op_len, op in aln.cigar if op in [0, 2, 3, 7, 8])
        
        assert q_consumed == 44, f"CIGAR should consume 44 query bases, got {q_consumed}"
        assert r_consumed == 43, f"CIGAR should consume 43 reference bases, got {r_consumed}"


def test_cigar_with_deletion():
    """Test that CIGAR correctly represents a 1-base deletion in the query."""
    # Reference has an extra 'A' at position 4 that's not in the query
    # Visual alignment:
    #   Reference: GAATAAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG
    #   Query:     GAAT-AGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG
    #              ^^^^D^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    reference = 'GAATAAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG'
    query = 'GAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG'
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write reference FASTA
        ref_file = os.path.join(tmpdir, 'ref.fa')
        with open(ref_file, 'w') as f:
            f.write('>ref\n')
            f.write(reference + '\n')
        
        # Build index
        indexer = BwaIndexer(verbose=0)
        index_prefix = indexer.build_index(ref_file)
        
        # Align query
        aligner = BwaAligner(index_prefix)
        alignments = aligner.align(query)
        
        # Should find one alignment
        assert len(alignments) == 1, f"Expected 1 alignment, got {len(alignments)}"
        
        aln = alignments[0]
        
        # Check basic alignment properties
        assert aln.q_st == 0
        assert aln.q_en == 43  # Full query aligned
        
        # Reference span should be 44 (longer than query by 1 due to deletion)
        ref_span = aln.r_en - aln.r_st
        assert ref_span == 44, f"Expected reference span of 44, got {ref_span}"
        
        # Check CIGAR string contains deletion operation
        assert 'D' in aln.cigar_str, f"CIGAR should contain deletion: {aln.cigar_str}"
        
        # Check specific CIGAR structure: should be 4M1D39M
        assert aln.cigar_str == '4M1D39M', f"Expected CIGAR '4M1D39M', got '{aln.cigar_str}'"
        
        # Check edit distance is 1 (one deletion)
        assert aln.NM == 1, f"Expected edit distance 1, got {aln.NM}"
        
        # Validate CIGAR consistency
        q_consumed = sum(op_len for op_len, op in aln.cigar if op in [0, 1, 4, 7, 8])
        r_consumed = sum(op_len for op_len, op in aln.cigar if op in [0, 2, 3, 7, 8])
        
        assert q_consumed == 43, f"CIGAR should consume 43 query bases, got {q_consumed}"
        assert r_consumed == 44, f"CIGAR should consume 44 reference bases, got {r_consumed}"


def test_cigar_perfect_match():
    """Test that CIGAR correctly represents a perfect match with no indels."""
    # Both sequences are identical - no insertions, deletions, or mismatches
    reference = 'GAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG'
    query = 'GAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG'
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write reference FASTA
        ref_file = os.path.join(tmpdir, 'ref.fa')
        with open(ref_file, 'w') as f:
            f.write('>ref\n')
            f.write(reference + '\n')
        
        # Build index
        indexer = BwaIndexer(verbose=0)
        index_prefix = indexer.build_index(ref_file)
        
        # Align query
        aligner = BwaAligner(index_prefix)
        alignments = aligner.align(query)
        
        # Should find one alignment
        assert len(alignments) == 1
        
        aln = alignments[0]
        
        # Check that it's a perfect match
        assert aln.q_st == 0
        assert aln.q_en == 43
        assert aln.r_en - aln.r_st == 43
        
        # CIGAR should be all matches
        assert aln.cigar_str == '43M', f"Expected CIGAR '43M', got '{aln.cigar_str}'"
        
        # Edit distance should be 0
        assert aln.NM == 0, f"Expected edit distance 0, got {aln.NM}"
