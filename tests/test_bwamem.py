"""Basic tests for bwamem package."""

import pytest
import bwamem


def test_bwa_aligner_import():
    """Test that BwaAligner can be imported."""
    assert hasattr(bwamem, "BwaAligner")
    assert callable(bwamem.BwaAligner)


def test_bwa_aligner_initialization():
    """Test BwaAligner initialization with invalid index."""
    # This should raise an error since we don't have a valid BWA index
    with pytest.raises(ValueError, match="Failed to load bwa index"):
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
    with pytest.raises(ValueError, match="Unknown algorithm"):
        bwamem.BwaIndexer(algorithm="invalid")


def test_bwa_indexer_algorithm_constants():
    """Test BwaIndexer algorithm constants."""
    assert bwamem.BwaIndexer.BWTALGO_AUTO == 0
    assert bwamem.BwaIndexer.BWTALGO_RB2 == 1
    assert bwamem.BwaIndexer.BWTALGO_BWTSW == 2
    assert bwamem.BwaIndexer.BWTALGO_IS == 3


def test_bwa_indexer_build_index_file_not_found():
    """Test BwaIndexer build_index with non-existent file."""
    indexer = bwamem.BwaIndexer()
    with pytest.raises(FileNotFoundError, match="FASTA file not found"):
        indexer.build_index("nonexistent.fasta")


def test_align_function_signature():
    """Test that the new align function has the correct signature."""
    from bwamem import BwaAligner, Alignment, PairedAlignment

    # Test that the classes exist
    assert Alignment is not None
    assert PairedAlignment is not None

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
    # No per-call insert parameters anymore


def test_data_structures():
    """Test that the data structures work correctly."""
    from bwamem import Alignment, PairedAlignment

    # Test Alignment
    se_aln = Alignment(
        ctg="chr1",
        ctg_len=10000,
        r_st=1000,
        strand=1,
        q_st=0,
        q_en=100,
        mapq=60,
        cigar=[[100, 0]],  # 100M as [length, op]
        NM=0,
        is_primary=True,
        read_num=0,
        trans_strand=0,
        score=100,
    )

    # Test core attributes
    assert se_aln.ctg == "chr1"
    assert se_aln.ctg_len == 10000
    assert se_aln.r_st == 1000
    assert se_aln.strand == 1
    assert se_aln.q_st == 0
    assert se_aln.q_en == 100
    assert se_aln.mapq == 60
    assert se_aln.cigar == [[100, 0]]
    assert se_aln.cigar_str == "100M"
    assert se_aln.NM == 0
    assert se_aln.score == 100
    assert se_aln.is_primary
    assert se_aln.read_num == 0
    assert se_aln.trans_strand == 0
    
    # Test calculated properties
    assert se_aln.r_en == 1100  # r_st + 100M
    assert se_aln.blen == 100
    assert se_aln.mlen == 100

    # Test PairedAlignment
    pe_aln = PairedAlignment(
        read1=se_aln, read2=se_aln, is_proper_pair=True, insert_size=500
    )

    assert pe_aln.read1 == se_aln
    assert pe_aln.read2 == se_aln
    assert pe_aln.is_proper_pair
    assert pe_aln.insert_size == 500


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
    # No per-call insert parameters anymore


def test_align_method_parameter_types():
    """Test that the align method accepts correct parameter types."""
    from bwamem import BwaAligner
    import inspect

    sig = inspect.signature(BwaAligner.align)

    # Check parameter types
    assert sig.parameters["seq1"].annotation is str
    assert sig.parameters["seq2"].annotation is str
    # Per-call insert parameters removed


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


def test_data_structure_fields():
    """Test that data structures have the correct fields."""
    from bwamem import Alignment, PairedAlignment

    # Test Alignment has the expected attributes
    test_aln = Alignment(
        ctg="test", ctg_len=1000, r_st=0, strand=1,
        q_st=0, q_en=100, mapq=60, cigar=[[100, 0]],
        NM=0, is_primary=True, read_num=0,
        trans_strand=0, score=100
    )
    
    # Test core attributes exist
    assert hasattr(test_aln, "ctg")
    assert hasattr(test_aln, "ctg_len")
    assert hasattr(test_aln, "r_st")
    assert hasattr(test_aln, "strand")
    assert hasattr(test_aln, "q_st")
    assert hasattr(test_aln, "q_en")
    assert hasattr(test_aln, "mapq")
    assert hasattr(test_aln, "cigar")
    assert hasattr(test_aln, "cigar_str")
    assert hasattr(test_aln, "NM")
    assert hasattr(test_aln, "is_primary")
    assert hasattr(test_aln, "read_num")
    assert hasattr(test_aln, "trans_strand")
    assert hasattr(test_aln, "score")
    
    # Test calculated properties
    assert hasattr(test_aln, "r_en")
    assert hasattr(test_aln, "blen")
    assert hasattr(test_aln, "mlen")

    # Test PairedAlignment fields
    expected_pe_fields = ("read1", "read2", "is_proper_pair", "insert_size")
    assert PairedAlignment._fields == expected_pe_fields
