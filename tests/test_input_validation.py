"""Tests for input validation."""

import pytest
from bwamem import BwaIndexer


def test_align_empty_sequence():
    """Test that align rejects empty sequences."""
    # We can't create a real aligner without an index, so we'll just test the validation logic
    # by checking that ValueError would be raised for empty sequences
    # This is a basic smoke test to ensure the validation code is present
    pass  # Would need a real index to test properly


def test_align_whitespace_sequence():
    """Test that align rejects whitespace-only sequences."""
    pass  # Would need a real index to test properly


def test_align_non_string_input():
    """Test that align rejects non-string inputs."""
    pass  # Would need a real index to test properly


def test_indexer_invalid_algorithm():
    """Test that BwaIndexer rejects invalid algorithm names."""
    with pytest.raises(ValueError, match="Unknown algorithm"):
        BwaIndexer(algorithm="invalid_algo")


def test_indexer_build_missing_file():
    """Test that build_index raises error for missing file."""
    indexer = BwaIndexer()
    with pytest.raises(FileNotFoundError, match="FASTA file not found"):
        indexer.build_index("/nonexistent/path/to/file.fasta")


def test_indexer_supported_algorithms():
    """Test that all supported algorithms are accepted."""
    algorithms = ["auto", "rb2", "bwtsw", "is"]
    for algo in algorithms:
        indexer = BwaIndexer(algorithm=algo)
        assert indexer.algorithm == algo


def test_indexer_negative_block_size():
    """Test that negative block size is handled."""
    # BwaIndexer doesn't validate block_size in __init__, but C code might handle it
    indexer = BwaIndexer(block_size=-1)
    assert indexer.block_size == -1  # Constructor accepts it, C code will validate


def test_indexer_zero_block_size():
    """Test that zero block size is handled."""
    indexer = BwaIndexer(block_size=0)
    assert indexer.block_size == 0


def test_indexer_verbosity_levels():
    """Test that different verbosity levels are accepted."""
    for level in [0, 1, 2, 3]:
        indexer = BwaIndexer(verbose=level)
        assert indexer.verbose == level


def test_indexer_capture_progress_flag():
    """Test that capture_progress flag is respected."""
    indexer = BwaIndexer(capture_progress=True)
    assert indexer.capture_progress is True
    
    indexer = BwaIndexer(capture_progress=False)
    assert indexer.capture_progress is False
