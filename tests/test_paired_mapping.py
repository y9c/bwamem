"""Test paired-end mapping logic."""

import os
import tempfile
import pytest
from bwamem import BwaIndexer, BwaAligner, read_paired_fastx, fastx_read

def test_paired_end_mapping():
    """Test paired-end mapping with synthetic data."""
    # Create synthetic reference and reads
    ref_seq = "ACGT" * 100
    r1_seq = "ACGT" * 5
    r2_seq = "ACGT" * 5
    
    with tempfile.TemporaryDirectory() as tmpdir:
        ref_path = os.path.join(tmpdir, "ref.fa")
        with open(ref_path, "w") as f:
            f.write(">ref\n" + ref_seq + "\n")
            
        indexer = BwaIndexer(verbose=0)
        idx_prefix = indexer.build_index(ref_path)
        aligner = BwaAligner(idx_prefix)
        
        # Test PE alignment
        results = aligner.align(r1_seq, r2_seq)
        assert isinstance(results, list)
        if results:
            for h1, h2, is_p, isize in results:
                if h1: assert isinstance(h1, tuple)
                if h2: assert isinstance(h2, tuple)

def test_single_end_mapping():
    """Test single-end mapping within the paired-end test file."""
    ref_seq = "ACGT" * 100
    r1_seq = "ACGT" * 5
    
    with tempfile.TemporaryDirectory() as tmpdir:
        ref_path = os.path.join(tmpdir, "ref.fa")
        with open(ref_path, "w") as f:
            f.write(">ref\n" + ref_seq + "\n")
            
        indexer = BwaIndexer(verbose=0)
        idx_prefix = indexer.build_index(ref_path)
        aligner = BwaAligner(idx_prefix)
        
        results = aligner.align(r1_seq)
        assert isinstance(results, list)
        if results:
            assert isinstance(results[0], tuple)
