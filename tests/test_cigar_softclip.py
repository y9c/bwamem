"""Test for correct CIGAR soft-clip operation (S instead of N).

This test validates that BWA correctly generates CIGAR strings with soft-clip
operations (S) instead of skip operations (N) when dealing with clipped reads.

Bug: Previously, the code in bwa/bwamem.c used operation code 3 (N - skip/intron)
instead of operation code 4 (S - soft-clip) when adding clipping to CIGAR strings.
This resulted in incorrect CIGARs like "16S16N27M" instead of "16S27M".

According to the BAM specification:
- M = 0 (match)
- I = 1 (insertion)
- D = 2 (deletion)
- N = 3 (skip/intron)
- S = 4 (soft-clip)
- H = 5 (hard-clip)
"""

import pytest
import tempfile
import os
from bwamem import BwaIndexer, BwaAligner


def test_cigar_softclip_no_skip_5prime():
    """Test that 5' soft-clips do not have adjacent N operations."""
    # Create a reference sequence
    reference = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'
    
    # Create a query with the last part matching but first part not matching
    # This should result in soft-clipping at the 5' end
    query = 'GGGGGGGGGGGGGGGGACGTACGTACGTACGTACGTACGT'
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write reference FASTA
        ref_file = os.path.join(tmpdir, 'ref.fa')
        with open(ref_file, 'w') as f:
            f.write('>chr1\n')
            f.write(reference + '\n')
        
        # Build index
        indexer = BwaIndexer(verbose=0)
        index_prefix = indexer.build_index(ref_file)
        
        # Align query
        aligner = BwaAligner(index_prefix)
        alignments = aligner.align(query)
        
        # Should find at least one alignment
        assert len(alignments) > 0, f"Expected at least 1 alignment, got {len(alignments)}"
        
        aln = alignments[0]
        
        # Check CIGAR string does NOT contain 'N' operation
        assert 'N' not in aln.cigar_str, f"CIGAR should not contain N (skip) operation: {aln.cigar_str}"
        
        # If there's soft-clipping at the start, verify the format
        if aln.cigar_str.startswith(('1S', '2S', '3S', '4S', '5S', '6S', '7S', '8S', '9S')) or 'S' in aln.cigar_str[:3]:
            # The CIGAR should be like "XS...M" without any N between S and M
            parts = aln.cigar_str.split('S')
            if len(parts) > 1:
                # After the S, there should be no N before the next operation
                remaining = parts[1]
                assert not remaining.startswith('N'), \
                    f"CIGAR should not have N after S: {aln.cigar_str}"


def test_cigar_softclip_no_skip_3prime():
    """Test that 3' soft-clips do not have adjacent N operations."""
    # Create a reference sequence
    reference = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'
    
    # Create a query with the first part matching but last part not matching
    # This should result in soft-clipping at the 3' end
    query = 'ACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGG'
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write reference FASTA
        ref_file = os.path.join(tmpdir, 'ref.fa')
        with open(ref_file, 'w') as f:
            f.write('>chr1\n')
            f.write(reference + '\n')
        
        # Build index
        indexer = BwaIndexer(verbose=0)
        index_prefix = indexer.build_index(ref_file)
        
        # Align query
        aligner = BwaAligner(index_prefix)
        alignments = aligner.align(query)
        
        # Should find at least one alignment
        assert len(alignments) > 0, f"Expected at least 1 alignment, got {len(alignments)}"
        
        aln = alignments[0]
        
        # Check CIGAR string does NOT contain 'N' operation
        assert 'N' not in aln.cigar_str, f"CIGAR should not contain N (skip) operation: {aln.cigar_str}"
        
        # If there's soft-clipping at the end, verify the format
        if 'S' in aln.cigar_str and aln.cigar_str.endswith('S'):
            # The CIGAR should be like "...MXS" without any N between M and S
            # Find the position before the last S
            parts = aln.cigar_str.rsplit('S', 1)
            if len(parts) > 1:
                before_s = parts[0]
                # Should not end with N
                assert not before_s.endswith('N'), \
                    f"CIGAR should not have N before final S: {aln.cigar_str}"


def test_cigar_softclip_both_ends():
    """Test that soft-clips on both ends do not have adjacent N operations."""
    # Create a reference sequence
    reference = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'
    
    # Create a query with mismatches at both ends
    # This should result in soft-clipping at both ends
    query = 'GGGGACGTACGTACGTACGTACGTACGTGGGG'
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write reference FASTA
        ref_file = os.path.join(tmpdir, 'ref.fa')
        with open(ref_file, 'w') as f:
            f.write('>chr1\n')
            f.write(reference + '\n')
        
        # Build index
        indexer = BwaIndexer(verbose=0)
        index_prefix = indexer.build_index(ref_file)
        
        # Align query
        aligner = BwaAligner(index_prefix)
        alignments = aligner.align(query)
        
        # Should find at least one alignment
        assert len(alignments) > 0, f"Expected at least 1 alignment, got {len(alignments)}"
        
        aln = alignments[0]
        
        # Check CIGAR string does NOT contain 'N' operation
        assert 'N' not in aln.cigar_str, f"CIGAR should not contain N (skip) operation: {aln.cigar_str}"


def test_cigar_operations_direct():
    """Test CIGAR operations directly by checking the cigar list structure."""
    # Create a simple test case
    reference = 'ACGTACGTACGTACGTACGTACGTACGTACGT'
    query = 'GGGGACGTACGTACGTACGTACGTGGGG'
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write reference FASTA
        ref_file = os.path.join(tmpdir, 'ref.fa')
        with open(ref_file, 'w') as f:
            f.write('>chr1\n')
            f.write(reference + '\n')
        
        # Build index
        indexer = BwaIndexer(verbose=0)
        index_prefix = indexer.build_index(ref_file)
        
        # Align query
        aligner = BwaAligner(index_prefix)
        alignments = aligner.align(query)
        
        # Should find at least one alignment
        assert len(alignments) > 0, f"Expected at least 1 alignment, got {len(alignments)}"
        
        aln = alignments[0]
        
        # Check the cigar list directly
        # CIGAR operations: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H
        for i, (length, op) in enumerate(aln.cigar):
            if op == 4:  # Soft-clip operation
                # Check that adjacent operations are not N (op=3)
                if i > 0:
                    prev_op = aln.cigar[i-1][1]
                    assert prev_op != 3, \
                        f"N operation (3) should not be adjacent to S operation (4): {aln.cigar}"
                if i < len(aln.cigar) - 1:
                    next_op = aln.cigar[i+1][1]
                    assert next_op != 3, \
                        f"N operation (3) should not be adjacent to S operation (4): {aln.cigar}"
            
            # Also verify that N operations don't appear incorrectly
            if op == 3:  # N (skip) operation
                # N operations should not be at the edges next to S
                if i > 0:
                    prev_op = aln.cigar[i-1][1]
                    assert prev_op != 4, \
                        f"N operation should not follow S operation: {aln.cigar}"
                if i < len(aln.cigar) - 1:
                    next_op = aln.cigar[i+1][1]
                    assert next_op != 4, \
                        f"N operation should not precede S operation: {aln.cigar}"
