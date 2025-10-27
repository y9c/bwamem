#!/usr/bin/env python3
"""
High-performance FASTQ reader optimized for speed.

This module provides fast reading of FASTQ files using optimized Python
with minimal overhead, similar to BWA's kseq approach.
"""

import gzip
import os
from typing import Iterator, Tuple, Optional
from collections import namedtuple

# Use namedtuple for better performance
FastqRead = namedtuple('FastqRead', ['name', 'comment', 'sequence', 'quality', 'length'])


class FastqReader:
    """
    High-performance FASTQ reader.
    
    This class provides fast reading of FASTQ files, both single-end and paired-end,
    using optimized Python implementation.
    """
    
    def __init__(self, file1: str, file2: Optional[str] = None):
        """
        Initialize FASTQ reader.
        
        Args:
            file1: Path to first FASTQ file (R1 for paired-end)
            file2: Path to second FASTQ file (R2 for paired-end, None for single-end)
        """
        self.file1 = file1
        self.file2 = file2
        self.is_paired = file2 is not None
        
        # Open files with appropriate compression handling
        self.fp1 = self._open_file(file1)
        self.fp2 = None
        if self.is_paired:
            self.fp2 = self._open_file(file2)
        
        self._closed = False
    
    def _open_file(self, file_path: str):
        """Open file with appropriate compression handling."""
        if file_path.endswith('.gz'):
            return gzip.open(file_path, 'rt')
        else:
            return open(file_path, 'r')
    
    def read(self) -> Union[FastqRead, Tuple[FastqRead, FastqRead], None]:
        """
        Read next sequence(s) from FASTQ file(s).
        
        Returns:
            For single-end: FastqRead object or None if end of file
            For paired-end: Tuple of (FastqRead, FastqRead) or None if end of file
        """
        if self._closed:
            return None
        
        try:
            if self.is_paired:
                # Read paired sequences
                line1 = self.fp1.readline()
                line2 = self.fp2.readline()
                
                if not line1 or not line2:
                    return None
                
                # Parse R1
                name1 = line1.strip()[1:]  # Remove '@'
                seq1 = self.fp1.readline().strip()
                self.fp1.readline()  # Skip '+'
                qual1 = self.fp1.readline().strip()
                
                # Parse R2
                name2 = line2.strip()[1:]  # Remove '@'
                seq2 = self.fp2.readline().strip()
                self.fp2.readline()  # Skip '+'
                qual2 = self.fp2.readline().strip()
                
                read1 = FastqRead(name1, "", seq1, qual1, len(seq1))
                read2 = FastqRead(name2, "", seq2, qual2, len(seq2))
                
                return (read1, read2)
            else:
                # Read single sequence
                line = self.fp1.readline()
                if not line:
                    return None
                
                name = line.strip()[1:]  # Remove '@'
                seq = self.fp1.readline().strip()
                self.fp1.readline()  # Skip '+'
                qual = self.fp1.readline().strip()
                
                return FastqRead(name, "", seq, qual, len(seq))
                
        except Exception:
            return None
    
    def close(self):
        """Close the FASTQ files."""
        if not self._closed:
            if self.fp1:
                self.fp1.close()
            if self.fp2:
                self.fp2.close()
            self._closed = True
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()
    
    def __iter__(self):
        """Make the reader iterable."""
        return self
    
    def __next__(self):
        """Get next sequence(s)."""
        result = self.read()
        if result is None:
            raise StopIteration
        return result


def read_fastq(file_path: str) -> Iterator[FastqRead]:
    """
    Read single-end FASTQ file.
    
    Args:
        file_path: Path to FASTQ file (supports .gz compression)
        
    Yields:
        FastqRead objects
    """
    with FastqReader(file_path) as reader:
        while True:
            read = reader.read()
            if read is None:
                break
            yield read


def read_paired_fastq(file1: str, file2: str) -> Iterator[Tuple[FastqRead, FastqRead]]:
    """
    Read paired-end FASTQ files.
    
    Args:
        file1: Path to R1 FASTQ file
        file2: Path to R2 FASTQ file
        
    Yields:
        Tuples of (FastqRead, FastqRead) objects
    """
    with FastqReader(file1, file2) as reader:
        while True:
            pair = reader.read()
            if pair is None:
                break
            yield pair


# Legacy compatibility - create a simple class for backward compatibility
class FastqReadLegacy:
    """Legacy FastqRead class for backward compatibility."""
    
    def __init__(self, qname: str, seq: str, qual: str):
        self.qname = qname
        self.seq = seq
        self.qual = qual
        self.name = qname  # Alias for compatibility
        self.sequence = seq  # Alias for compatibility
        self.quality = qual  # Alias for compatibility
        self.length = len(seq)  # Alias for compatibility
    
    def __repr__(self):
        return f"FastqRead(qname='{self.qname}', seq_len={len(self.seq)}, qual_len={len(self.qual)})"