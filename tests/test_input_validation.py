"""Tests for input validation and error handling."""

import pytest
import os
import tempfile
from bwamem import BwaIndexer, BwaAligner


def test_align_empty_sequence():
    """Test alignment with empty sequence."""
    # This shouldn't crash, but return no hits
    # Need a valid aligner first
    with tempfile.NamedTemporaryFile(suffix=".fa", mode="w", delete=False) as f:
        f.write(">ref\nACGTACGTACGTACGTACGTACGTACGTACGT\n")
        ref_path = f.name
    
    try:
        indexer = BwaIndexer(verbose=0)
        idx = indexer.build_index(ref_path)
        aligner = BwaAligner(idx)
        
        assert aligner.align("") == []
        assert aligner.align("", "ACGT") == []
    finally:
        if os.path.exists(ref_path): os.remove(ref_path)


def test_align_whitespace_sequence():
    """Test alignment with whitespace sequence."""
    with tempfile.NamedTemporaryFile(suffix=".fa", mode="w", delete=False) as f:
        f.write(">ref\nACGTACGTACGTACGTACGTACGTACGTACGT\n")
        ref_path = f.name
    
    try:
        indexer = BwaIndexer(verbose=0)
        idx = indexer.build_index(ref_path)
        aligner = BwaAligner(idx)
        
        assert aligner.align("   ") == []
    finally:
        if os.path.exists(ref_path): os.remove(ref_path)


def test_align_non_string_input():
    """Test alignment with non-string input."""
    # This should handle the error gracefully or raise TypeError
    with tempfile.NamedTemporaryFile(suffix=".fa", mode="w", delete=False) as f:
        f.write(">ref\nACGTACGTACGTACGTACGTACGTACGTACGT\n")
        ref_path = f.name
    
    try:
        indexer = BwaIndexer(verbose=0)
        idx = indexer.build_index(ref_path)
        aligner = BwaAligner(idx)
        
        with pytest.raises((TypeError, AttributeError)):
            aligner.align(123)
    finally:
        if os.path.exists(ref_path): os.remove(ref_path)


def test_indexer_invalid_algorithm():
    """Test indexer with invalid algorithm."""
    with pytest.raises(KeyError):
        BwaIndexer(algorithm="invalid_algo")


def test_indexer_build_missing_file():
    """Test building index with missing file."""
    indexer = BwaIndexer()
    with pytest.raises(FileNotFoundError):
        indexer.build_index("/nonexistent/path/to/file.fasta")


def test_indexer_supported_algorithms():
    """Test that all supported algorithms work."""
    for algo in ["auto", "rb2", "bwtsw", "is"]:
        indexer = BwaIndexer(algorithm=algo)
        assert indexer.algorithm == algo


def test_indexer_negative_block_size():
    """Test indexer with negative block size."""
    # Block size is used in C, should be positive
    indexer = BwaIndexer(block_size=-1)
    assert indexer.block_size == -1


def test_indexer_zero_block_size():
    """Test indexer with zero block size."""
    indexer = BwaIndexer(block_size=0)
    assert indexer.block_size == 0


def test_indexer_verbosity_levels():
    """Test different verbosity levels."""
    for level in [0, 1, 2, 3]:
        indexer = BwaIndexer(verbose=level)
        assert indexer.verbose == level


def test_indexer_capture_progress_flag():
    """Test capture_progress flag."""
    indexer = BwaIndexer(capture_progress=True)
    assert indexer.capture_progress is True
    
    indexer = BwaIndexer(capture_progress=False)
    assert indexer.capture_progress is False
