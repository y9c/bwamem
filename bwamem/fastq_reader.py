#!/usr/bin/env python3
"""
Helper functions for reading FASTQ files sequentially.
Based on the BWA kseq implementation.
"""

import gzip
from typing import Iterator, Tuple, Optional, TextIO


class FastqRead:
    """Represents a single FASTQ read with qname, sequence, and quality."""
    
    def __init__(self, qname: str, seq: str, qual: str):
        self.qname = qname
        self.seq = seq
        self.qual = qual
    
    def __repr__(self):
        return f"FastqRead(qname='{self.qname}', seq_len={len(self.seq)}, qual_len={len(self.qual)})"


def read_fastq(file_path: str) -> Iterator[FastqRead]:
    """
    Read FASTQ file sequentially and yield FastqRead objects.
    
    Args:
        file_path: Path to FASTQ file (supports .gz compression)
    
    Yields:
        FastqRead objects containing qname, seq, and qual
    """
    # Open file with gzip support
    if file_path.endswith('.gz'):
        file_handle = gzip.open(file_path, 'rt')
    else:
        file_handle = open(file_path, 'r')
    
    try:
        while True:
            # Read header line (starts with @)
            header = file_handle.readline().strip()
            if not header:  # End of file
                break
            if not header.startswith('@'):
                raise ValueError(f"Invalid FASTQ format: expected '@' at start of header, got '{header[:10]}...'")
            
            # Extract qname (remove @ and any /1 or /2 suffix)
            qname = header[1:].rstrip('/12')
            
            # Read sequence line
            seq = file_handle.readline().strip()
            if not seq:
                raise ValueError("Unexpected end of file while reading sequence")
            
            # Read + line (should start with +)
            plus_line = file_handle.readline().strip()
            if not plus_line.startswith('+'):
                raise ValueError(f"Invalid FASTQ format: expected '+' line, got '{plus_line[:10]}...'")
            
            # Read quality line
            qual = file_handle.readline().strip()
            if not qual:
                raise ValueError("Unexpected end of file while reading quality")
            
            # Validate sequence and quality lengths match
            if len(seq) != len(qual):
                raise ValueError(f"Sequence and quality lengths don't match: {len(seq)} vs {len(qual)}")
            
            yield FastqRead(qname, seq, qual)
    
    finally:
        file_handle.close()


def read_paired_fastq(file1_path: str, file2_path: str) -> Iterator[Tuple[FastqRead, FastqRead]]:
    """
    Read paired FASTQ files and yield pairs of FastqRead objects.
    
    Args:
        file1_path: Path to first FASTQ file (R1)
        file2_path: Path to second FASTQ file (R2)
    
    Yields:
        Tuples of (FastqRead, FastqRead) for paired reads
    """
    read1_iter = read_fastq(file1_path)
    read2_iter = read_fastq(file2_path)
    
    try:
        while True:
            read1 = next(read1_iter, None)
            read2 = next(read2_iter, None)
            
            if read1 is None and read2 is None:
                break  # Both files exhausted
            elif read1 is None:
                raise ValueError("R1 file ended before R2 file")
            elif read2 is None:
                raise ValueError("R2 file ended before R1 file")
            
            # Validate that qnames match (should be identical for paired reads)
            if read1.qname != read2.qname:
                raise ValueError(f"Paired read qnames don't match: '{read1.qname}' vs '{read2.qname}'")
            
            yield (read1, read2)
    
    except StopIteration:
        pass
