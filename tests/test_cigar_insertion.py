"""Test for correct CIGAR generation with insertions.

This test validates that BWA correctly generates CIGAR strings for alignments
with insertions, ensuring that the reference span is calculated correctly.
"""

import pytest
import tempfile
import os
import re
from bwamem import BwaIndexer, BwaAligner


def parse_cigar(cigar_str):
    return re.findall(r"(\d+)([MIDNSHP=X])", cigar_str)


def test_cigar_with_insertion():
    """Test that CIGAR correctly represents a 1-base insertion in the query."""
    prefix = "GGCGAGCCACCGCCCGTCCCCGCCCCTTGCCTCTCGGCGCCCCCTCGATGCTCTTAGCTGAGTGTCCCGCGGGGCCCGAAGCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCCCGAGCCGCCTGGATACCGCAGCTAGGAATAATGGAAT"
    alignment_region = "GAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG"
    suffix = "CCATGATTAAGAGGGACGGCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACCAGAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTC"
    reference = prefix + alignment_region + suffix
    query = "GAATAAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG"

    with tempfile.TemporaryDirectory() as tmpdir:
        ref_file = os.path.join(tmpdir, "ref.fa")
        with open(ref_file, "w") as f:
            f.write(">ref\n" + reference + "\n")

        indexer = BwaIndexer(verbose=0)
        index_prefix = indexer.build_index(ref_file)

        aligner = BwaAligner(index_prefix)
        alignments = aligner.align(query)

        assert len(alignments) == 1
        aln = alignments[0]

        # (ctg, r_st, r_en, strand, q_st, q_en, mapq, cigar_str, NM, score)
        assert aln[0] == "ref"
        assert aln[4] == 0
        assert aln[5] == 44
        assert aln[1] == len(prefix)

        ref_span = aln[2] - aln[1]
        assert ref_span == 43

        assert "I" in aln[7]
        assert aln[7] == "4M1I39M"
        assert aln[8] == 1

        cigar_ops = parse_cigar(aln[7])
        q_consumed = sum(int(l) for l, op in cigar_ops if op in "MIS=X")
        r_consumed = sum(int(l) for l, op in cigar_ops if op in "MDN=X")

        assert q_consumed == 44
        assert r_consumed == 43


def test_cigar_with_deletion():
    """Test that CIGAR correctly represents a 1-base deletion in the query."""
    reference = "GAATAAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG"
    query = "GAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG"

    with tempfile.TemporaryDirectory() as tmpdir:
        ref_file = os.path.join(tmpdir, "ref.fa")
        with open(ref_file, "w") as f:
            f.write(">ref\n" + reference + "\n")

        indexer = BwaIndexer(verbose=0)
        index_prefix = indexer.build_index(ref_file)

        aligner = BwaAligner(index_prefix)
        alignments = aligner.align(query)

        assert len(alignments) == 1
        aln = alignments[0]

        assert aln[4] == 0
        assert aln[5] == 43

        ref_span = aln[2] - aln[1]
        assert ref_span == 44

        assert "D" in aln[7]
        assert aln[7] == "4M1D39M"
        assert aln[8] == 1

        cigar_ops = parse_cigar(aln[7])
        q_consumed = sum(int(l) for l, op in cigar_ops if op in "MIS=X")
        r_consumed = sum(int(l) for l, op in cigar_ops if op in "MDN=X")

        assert q_consumed == 43
        assert r_consumed == 44


def test_cigar_perfect_match():
    """Test that CIGAR correctly represents a perfect match with no indels."""
    reference = "GAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG"
    query = "GAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG"

    with tempfile.TemporaryDirectory() as tmpdir:
        ref_file = os.path.join(tmpdir, "ref.fa")
        with open(ref_file, "w") as f:
            f.write(">ref\n" + reference + "\n")

        indexer = BwaIndexer(verbose=0)
        index_prefix = indexer.build_index(ref_file)

        aligner = BwaAligner(index_prefix)
        alignments = aligner.align(query)

        assert len(alignments) == 1
        aln = alignments[0]

        assert aln[4] == 0
        assert aln[5] == 43
        assert aln[2] - aln[1] == 43
        assert aln[7] == "43M"
        assert aln[8] == 0
