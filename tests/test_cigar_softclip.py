"""Test for correct CIGAR generation with soft-clipping.

This test validates that BWA correctly generates CIGAR strings with soft-clipping
when the alignment doesn't start or end at the query boundaries.
"""

import pytest
import tempfile
import os
import re
from bwamem import BwaIndexer, BwaAligner

def parse_cigar(cigar_str):
    return re.findall(r'(\d+)([MIDNSHP=X])', cigar_str)

def test_cigar_softclip_no_skip_5prime():
    """Test that 5' soft-clipping is correctly represented."""
    reference = 'GAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG'
    # Query has an extra 'CCCCC' at the beginning
    query = 'CCCCCGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG'
    
    with tempfile.TemporaryDirectory() as tmpdir:
        ref_file = os.path.join(tmpdir, 'ref.fa')
        with open(ref_file, 'w') as f:
            f.write('>ref\n' + reference + '\n')
        
        indexer = BwaIndexer(verbose=0)
        index_prefix = indexer.build_index(ref_file)
        
        aligner = BwaAligner(index_prefix)
        alignments = aligner.align(query)
        
        assert len(alignments) == 1
        aln = alignments[0]
        
        # (ctg, r_st, r_en, strand, q_st, q_en, mapq, cigar_str, NM, score)
        assert 'S' in aln[7]
        assert 'N' not in aln[7]
        
        # Validate total query length
        cigar_ops = parse_cigar(aln[7])
        q_len = sum(int(l) for l, op in cigar_ops if op in 'MIS=X')
        assert q_len == len(query)


def test_cigar_softclip_no_skip_3prime():
    """Test that 3' soft-clipping is correctly represented."""
    reference = 'GAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG'
    # Query has an extra 'GGGGG' at the end
    query = 'GAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAAAAAGGGGG'
    
    with tempfile.TemporaryDirectory() as tmpdir:
        ref_file = os.path.join(tmpdir, 'ref.fa')
        with open(ref_file, 'w') as f:
            f.write('>ref\n' + reference + '\n')
        
        indexer = BwaIndexer(verbose=0)
        index_prefix = indexer.build_index(ref_file)
        
        aligner = BwaAligner(index_prefix)
        alignments = aligner.align(query)
        
        assert len(alignments) == 1
        aln = alignments[0]
        
        assert 'S' in aln[7]
        assert 'N' not in aln[7]
        
        cigar_ops = parse_cigar(aln[7])
        q_len = sum(int(l) for l, op in cigar_ops if op in 'MIS=X')
        assert q_len == len(query)


def test_cigar_softclip_both_ends():
    """Test that soft-clipping on both ends is correctly represented."""
    reference = 'GAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG'
    query = 'AAAAAGAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACCCCCC'
    
    with tempfile.TemporaryDirectory() as tmpdir:
        ref_file = os.path.join(tmpdir, 'ref.fa')
        with open(ref_file, 'w') as f:
            f.write('>ref\n' + reference + '\n')
        
        indexer = BwaIndexer(verbose=0)
        index_prefix = indexer.build_index(ref_file)
        
        aligner = BwaAligner(index_prefix)
        alignments = aligner.align(query)
        
        assert len(alignments) == 1
        aln = alignments[0]
        
        assert 'S' in aln[7]
        assert 'N' not in aln[7]
        
        cigar_ops = parse_cigar(aln[7])
        q_len = sum(int(l) for l, op in cigar_ops if op in 'MIS=X')
        assert q_len == len(query)


def test_cigar_operations_direct():
    """Test direct access to CIGAR operations."""
    reference = 'GAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGG'
    query = 'AAAAAGAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACCCCCC'
    
    with tempfile.TemporaryDirectory() as tmpdir:
        ref_file = os.path.join(tmpdir, 'ref.fa')
        with open(ref_file, 'w') as f:
            f.write('>ref\n' + reference + '\n')
        
        indexer = BwaIndexer(verbose=0)
        index_prefix = indexer.build_index(ref_file)
        
        aligner = BwaAligner(index_prefix)
        alignments = aligner.align(query)
        
        aln = alignments[0]
        cigar_ops = parse_cigar(aln[7])
        
        # Verify it starts/ends with S
        assert cigar_ops[0][1] == 'S'
        assert cigar_ops[-1][1] == 'S'
        
        q_len = sum(int(l) for l, op in cigar_ops if op in 'MIS=X')
        assert q_len == len(query)
