"""
Tests for FASTA/FASTQ reader functionality.
"""

import os
import tempfile
import gzip
from bwamem import fastx_read, FastxReader


def create_test_fasta(path):
    """Create a test FASTA file."""
    with open(path, "w") as f:
        f.write(">seq1 test sequence 1\n")
        f.write("ACGTACGTACGTACGT\n")
        f.write("ACGTACGTACGTACGT\n")
        f.write(">seq2 test sequence 2\n")
        f.write("TGCATGCATGCATGCA\n")
        f.write(">seq3\n")
        f.write("GGGGCCCCAAAATTTT\n")


def create_test_fastq(path):
    """Create a test FASTQ file."""
    with open(path, "w") as f:
        f.write("@seq1 test sequence 1\n")
        f.write("ACGTACGTACGTACGT\n")
        f.write("+\n")
        f.write("IIIIIIIIIIIIIIII\n")
        f.write("@seq2 test sequence 2\n")
        f.write("TGCATGCATGCATGCA\n")
        f.write("+\n")
        f.write("HHHHHHHHHHHHHHHH\n")
        f.write("@seq3\n")
        f.write("GGGGCCCCAAAATTTT\n")
        f.write("+\n")
        f.write("FFFFFFFFFFFFFFFF\n")


def test_read_fasta_format():
    """Test that fastx_read can read FASTA format."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        fasta_path = f.name
        create_test_fasta(fasta_path)

    try:
        reads = list(fastx_read(fasta_path))

        # Check number of reads
        assert len(reads) == 3, f"Expected 3 reads, got {len(reads)}"

        # Check first read
        assert reads[0].name == "seq1"
        assert reads[0].comment == "test sequence 1"
        assert reads[0].sequence == "ACGTACGTACGTACGTACGTACGTACGTACGT"
        assert reads[0].quality == "", "FASTA should have empty quality string"
        assert reads[0].length == 32

        # Check second read
        assert reads[1].name == "seq2"
        assert reads[1].sequence == "TGCATGCATGCATGCA"
        assert reads[1].quality == "", "FASTA should have empty quality string"

        # Check all reads have empty quality
        assert all(len(r.quality) == 0 for r in reads), (
            "All FASTA reads should have empty quality"
        )

    finally:
        os.unlink(fasta_path)


def test_read_fastq_format():
    """Test that fastx_read can read FASTQ format."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fastq", delete=False) as f:
        fastq_path = f.name
        create_test_fastq(fastq_path)

    try:
        reads = list(fastx_read(fastq_path))

        # Check number of reads
        assert len(reads) == 3, f"Expected 3 reads, got {len(reads)}"

        # Check first read
        assert reads[0].name == "seq1"
        assert reads[0].comment == "test sequence 1"
        assert reads[0].sequence == "ACGTACGTACGTACGT"
        assert reads[0].quality == "IIIIIIIIIIIIIIII", (
            "FASTQ should have quality string"
        )
        assert reads[0].length == 16

        # Check all reads have quality strings
        assert all(len(r.quality) > 0 for r in reads), (
            "All FASTQ reads should have quality strings"
        )

    finally:
        os.unlink(fastq_path)


def test_read_gzipped_fasta():
    """Test that fastx_read can read gzip-compressed FASTA files."""
    with tempfile.NamedTemporaryFile(suffix=".fasta.gz", delete=False) as f:
        gz_path = f.name

    # Create gzipped FASTA
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        fasta_path = f.name
        create_test_fasta(fasta_path)

    try:
        # Gzip the file
        with open(fasta_path, "rb") as f_in:
            with gzip.open(gz_path, "wb") as f_out:
                f_out.write(f_in.read())

        reads = list(fastx_read(gz_path))

        assert len(reads) == 3
        assert all(len(r.quality) == 0 for r in reads), (
            "Gzipped FASTA should have empty quality"
        )

    finally:
        os.unlink(fasta_path)
        os.unlink(gz_path)


def test_read_gzipped_fastq():
    """Test that fastx_read can read gzip-compressed FASTQ files."""
    with tempfile.NamedTemporaryFile(suffix=".fastq.gz", delete=False) as f:
        gz_path = f.name

    # Create gzipped FASTQ
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fastq", delete=False) as f:
        fastq_path = f.name
        create_test_fastq(fastq_path)

    try:
        # Gzip the file
        with open(fastq_path, "rb") as f_in:
            with gzip.open(gz_path, "wb") as f_out:
                f_out.write(f_in.read())

        reads = list(fastx_read(gz_path))

        assert len(reads) == 3
        assert all(len(r.quality) > 0 for r in reads), (
            "Gzipped FASTQ should have quality strings"
        )

    finally:
        os.unlink(fastq_path)
        os.unlink(gz_path)


def test_fastx_reader_context_manager():
    """Test that FastxReader works as a context manager."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        fasta_path = f.name
        create_test_fasta(fasta_path)

    try:
        count = 0
        with FastxReader(fasta_path) as reader:
            for read in reader:
                count += 1

        assert count == 3, f"Expected 3 reads, got {count}"

    finally:
        os.unlink(fasta_path)


def test_fastx_reader_iterator():
    """Test that FastxReader is iterable."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fastq", delete=False) as f:
        fastq_path = f.name
        create_test_fastq(fastq_path)

    try:
        reader = FastxReader(fastq_path)
        reads = list(reader)
        reader.close()

        assert len(reads) == 3
        assert all(hasattr(r, "sequence") for r in reads)

    finally:
        os.unlink(fastq_path)
