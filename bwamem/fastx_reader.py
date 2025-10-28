#!/usr/bin/env python3
"""
High-performance FASTA/FASTQ reader using BWA's kseq library directly.

This module provides the fastest possible FASTA/FASTQ reading by using BWA's
optimized C implementation directly through cffi. Supports both FASTA and FASTQ
formats with automatic format detection.
"""

from collections import namedtuple
from typing import Iterator, Optional, Tuple, Union

from .libbwa import ffi, libbwa

# Use namedtuple for better performance
FastxRead = namedtuple(
    "FastxRead", ["name", "comment", "sequence", "quality", "length"]
)


class FastxReader:
    """
    High-performance FASTA/FASTQ reader using BWA's kseq library directly.

    This class provides the fastest possible reading of FASTA/FASTQ files by using
    BWA's optimized C implementation directly. Supports both formats with automatic
    format detection.
    """

    def __init__(self, file1: str, file2: Optional[str] = None):
        """
        Initialize FASTA/FASTQ reader.

        Args:
            file1: Path to first FASTA/FASTQ file (R1 for paired-end)
            file2: Path to second FASTA/FASTQ file (R2 for paired-end, None for single-end)
        """
        self.file1 = file1
        self.file2 = file2
        self.is_paired = file2 is not None

        # Open files using BWA's err_xzopen_core (handles gzip automatically)
        self.fp1 = libbwa.err_xzopen_core(b"FastxReader", file1.encode("utf-8"), b"r")
        if not self.fp1:
            raise IOError(f"Failed to open file: {file1}")

        self.fp2 = None
        if self.is_paired:
            self.fp2 = libbwa.err_xzopen_core(
                b"FastxReader", file2.encode("utf-8"), b"r"
            )
            if not self.fp2:
                libbwa.gzclose(self.fp1)
                raise IOError(f"Failed to open file: {file2}")

        # Initialize kseq readers
        self.ks1 = libbwa.kseq_init(self.fp1)
        if not self.ks1:
            libbwa.gzclose(self.fp1)
            if self.fp2:
                libbwa.gzclose(self.fp2)
            raise IOError("Failed to initialize kseq for first file")

        self.ks2 = None
        if self.is_paired:
            self.ks2 = libbwa.kseq_init(self.fp2)
            if not self.ks2:
                libbwa.kseq_destroy(self.ks1)
                libbwa.gzclose(self.fp1)
                libbwa.gzclose(self.fp2)
                raise IOError("Failed to initialize kseq for second file")

        self._closed = False

    def read(self) -> Union[FastxRead, Tuple[FastxRead, FastxRead], None]:
        """
        Read next sequence(s) from FASTA/FASTQ file(s).

        Returns:
            For single-end: FastxRead object or None if end of file
            For paired-end: Tuple of (FastxRead, FastxRead) or None if end of file
        """
        if self._closed:
            return None

        try:
            if self.is_paired:
                # Use bseq_read for paired-end processing
                n = ffi.new("int *")
                seqs = libbwa.bseq_read(2, n, self.ks1, self.ks2)

                if not seqs or n[0] < 2:
                    if seqs:
                        libbwa.free(seqs)
                    return None

                # Paired-end: return first two sequences
                seq1 = self._bseq_to_fastx_read(seqs[0])
                seq2 = self._bseq_to_fastx_read(seqs[1])
                libbwa.free(seqs)
                return (seq1, seq2)
            else:
                # Single-end: use kseq_read directly
                ret = libbwa.kseq_read(self.ks1)
                if ret < 0:
                    return None

                # Convert kseq_t to FastxRead
                return self._kseq_to_fastx_read(self.ks1)

        except Exception:
            return None

    def _bseq_to_fastx_read(self, bseq):
        """Convert bseq1_t to FastxRead."""
        # Access bseq1_t fields
        name = ffi.string(bseq.name).decode("utf-8") if bseq.name else ""
        comment = ffi.string(bseq.comment).decode("utf-8") if bseq.comment else ""
        seq = ffi.string(bseq.seq).decode("utf-8") if bseq.seq else ""
        qual = ffi.string(bseq.qual).decode("utf-8") if bseq.qual else ""
        length = bseq.l_seq

        return FastxRead(name, comment, seq, qual, length)

    def _kseq_to_fastx_read(self, ks):
        """Convert kseq_t to FastxRead."""
        # Access kseq_t fields (kstring_t structures)
        name = ffi.string(ks.name.s).decode("utf-8") if ks.name.s else ""
        comment = ffi.string(ks.comment.s).decode("utf-8") if ks.comment.s else ""
        seq = ffi.string(ks.seq.s).decode("utf-8") if ks.seq.s else ""
        qual = ffi.string(ks.qual.s).decode("utf-8") if ks.qual.s else ""
        length = ks.seq.l

        return FastxRead(name, comment, seq, qual, length)

    def close(self):
        """Close the FASTA/FASTQ files."""
        if not self._closed:
            if self.ks2:
                libbwa.kseq_destroy(self.ks2)
                libbwa.gzclose(self.fp2)
            if self.ks1:
                libbwa.kseq_destroy(self.ks1)
                libbwa.gzclose(self.fp1)
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


def fastx_read(file_path: str) -> Iterator[FastxRead]:
    """
    Read FASTA or FASTQ file using BWA C implementation.

    Supports both FASTA and FASTQ formats, with automatic format detection.
    For FASTA files, the quality field will be empty.

    Args:
        file_path: Path to FASTA/FASTQ file (supports .gz compression)

    Yields:
        FastxRead objects
    """
    with FastxReader(file_path) as reader:
        while True:
            read = reader.read()
            if read is None:
                break
            yield read


def read_paired_fastx(file1: str, file2: str) -> Iterator[Tuple[FastxRead, FastxRead]]:
    """
    Read paired-end FASTA/FASTQ files using BWA C implementation.

    Args:
        file1: Path to R1 FASTA/FASTQ file
        file2: Path to R2 FASTA/FASTQ file

    Yields:
        Tuples of (FastxRead, FastxRead) objects
    """
    with FastxReader(file1, file2) as reader:
        while True:
            pair = reader.read()
            if pair is None:
                break
            yield pair
