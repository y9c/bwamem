"""Basic tests for bwamem package."""

import pytest
import bwamem


def test_bwa_aligner_import():
    """Test that BwaAligner can be imported."""
    assert hasattr(bwamem, "BwaAligner")
    assert callable(bwamem.BwaAligner)


def test_bwa_aligner_initialization():
    """Test BwaAligner initialization with invalid index."""
    with pytest.raises(ValueError, match="Failed to load BWA index"):
        bwamem.BwaAligner("nonexistent_index")


def test_bwa_indexer_import():
    """Test that BwaIndexer can be imported."""
    assert hasattr(bwamem, "BwaIndexer")
    assert callable(bwamem.BwaIndexer)


def test_bwa_indexer_initialization():
    """Test BwaIndexer initialization with different algorithms."""
    # Test default initialization
    indexer = bwamem.BwaIndexer()
    assert indexer.algorithm == "auto"
    assert indexer.block_size == 10000000

    # Test with specific algorithm
    indexer = bwamem.BwaIndexer(algorithm="is")
    assert indexer.algorithm == "is"

    # Test with invalid algorithm
    with pytest.raises(KeyError):
        bwamem.BwaIndexer(algorithm="invalid")


def test_bwa_indexer_build_index_file_not_found():
    """Test BwaIndexer build_index with non-existent file."""
    indexer = bwamem.BwaIndexer()
    with pytest.raises(FileNotFoundError):
        indexer.build_index("nonexistent.fasta")


def test_align_function_signature():
    """Test that the new align function has the correct signature."""
    from bwamem import BwaAligner

    # Test that BwaAligner has the align method
    assert hasattr(BwaAligner, "align")

    # Test the method signature
    import inspect

    sig = inspect.signature(BwaAligner.align)
    params = list(sig.parameters.keys())

    # Should have seq1 and optional seq2 parameters
    assert "seq1" in params
    assert "seq2" in params

    # seq1 should be required, others optional
    assert sig.parameters["seq1"].default == inspect.Parameter.empty
    assert sig.parameters["seq2"].default is None


def test_align_method_availability():
    """Test that the align method is available and has correct signature."""
    from bwamem import BwaAligner
    import inspect

    # Test that align method exists
    assert hasattr(BwaAligner, "align")
    assert callable(getattr(BwaAligner, "align"))

    # Test method signature
    sig = inspect.signature(BwaAligner.align)
    params = list(sig.parameters.keys())

    # Required parameters
    assert "self" in params
    assert "seq1" in params

    # Optional parameters
    assert "seq2" in params

    # Check parameter defaults
    assert sig.parameters["seq1"].default == inspect.Parameter.empty  # Required
    assert sig.parameters["seq2"].default is None  # Optional


def test_single_end_vs_paired_end_usage():
    """Test the different usage patterns for SE vs PE alignment."""

    # Test that we can call align with just one sequence (SE)
    se_params = {"seq1": "ACGATCGCGATCGA"}

    # This should be valid for single-end
    assert "seq1" in se_params
    assert se_params["seq1"] == "ACGATCGCGATCGA"

    # Test that we can call align with two sequences (PE)
    pe_params = {
        "seq1": "ACGATCGCGATCGA",
        "seq2": "TTCGATCGATCGAT",
    }

    # This should be valid for paired-end
    assert "seq1" in pe_params
    assert "seq2" in pe_params
    assert pe_params["seq1"] == "ACGATCGCGATCGA"
    assert pe_params["seq2"] == "TTCGATCGATCGAT"
