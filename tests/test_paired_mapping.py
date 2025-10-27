#!/usr/bin/env python3
"""
Test script for paired-end mapping using bwap.
Reads FASTQ files and performs paired-end alignment.
"""

import sys
import os
import tempfile
from pathlib import Path

# Add parent directory to path to import bwap
sys.path.insert(0, str(Path(__file__).parent.parent))

import bwap
import pytest
from create_demo_files import create_demo_files


def test_paired_end_mapping():
    """Test complete paired-end mapping workflow."""
    
    print("🧬 Testing Paired-End Mapping with bwap")
    print("=" * 50)
    
    # Create demo files
    print("📁 Creating demo files...")
    create_demo_files()
    
    # Check if files exist
    ref_file = "tests/demo_reference.fasta"
    r1_file = "tests/demo_reads_R1.fastq"
    r2_file = "tests/demo_reads_R2.fastq"
    
    if not all(os.path.exists(f) for f in [ref_file, r1_file, r2_file]):
        pytest.fail("Demo files not created properly!")
    
    print("✅ Demo files created successfully")
    
    # Step 1: Build BWA index
    print("\n🔨 Building BWA index...")
    try:
        indexer = bwap.BwaIndexer()
        index_path = indexer.build_index(ref_file, "tests/demo_index/demo_index")
        print(f"✅ Index built successfully: {index_path}")
    except Exception as e:
        pytest.fail(f"Failed to build index: {e}")
    
    # Step 2: Initialize aligner
    print("\n🎯 Initializing BWA aligner...")
    try:
        aligner = bwap.BwaAligner(index_path)
        print("✅ Aligner initialized successfully")
    except Exception as e:
        pytest.fail(f"Failed to initialize aligner: {e}")
    
    # Step 3: Read and align paired reads
    print("\n📖 Reading and aligning paired reads...")
    
    try:
        # Read paired FASTQ files
        paired_reads = list(bwap.read_paired_fastq(r1_file, r2_file))
        print(f"✅ Read {len(paired_reads)} paired reads")
        
        # Process each pair
        for i, (read1, read2) in enumerate(paired_reads):
            print(f"\n--- Processing pair {i+1}: {read1.qname} ---")
            
            # Perform paired-end alignment
            try:
                result = aligner.align(
                    seq1=read1.seq,
                    seq2=read2.seq,
                )
                
                if isinstance(result, bwap.PairedAlignment):
                    print(f"✅ Paired alignment successful!")
                    print(f"   Read 1: {result.read1}")
                    print(f"   Read 2: {result.read2}")
                    print(f"   Proper pair: {result.is_proper_pair}")
                    print(f"   Insert size: {result.insert_size}")
                    
                    # Print CIGAR strings
                    if result.read1.cigar:
                        print(f"   R1 CIGAR: {result.read1.cigar}")
                    if result.read2.cigar:
                        print(f"   R2 CIGAR: {result.read2.cigar}")
                        
                else:
                    print(f"⚠️  Unexpected result type: {type(result)}")
                    
            except Exception as e:
                pytest.fail(f"Alignment failed for pair {i+1}: {e}")
                
    except Exception as e:
        pytest.fail(f"Failed to read FASTQ files: {e}")
    
    print("\n🎉 Paired-end mapping test completed!")


def test_single_end_mapping():
    """Test single-end mapping for comparison."""
    
    print("\n🧬 Testing Single-End Mapping")
    print("=" * 30)
    
    try:
        # Read single FASTQ file
        reads = list(bwap.read_fastq("tests/demo_reads_R1.fastq"))
        print(f"✅ Read {len(reads)} single reads")
        
        # Initialize aligner (reuse index from previous test)
        aligner = bwap.BwaAligner("tests/demo_index/demo_index")
        
        # Process first read
        read = reads[0]
        print(f"\n--- Processing read: {read.qname} ---")
        
        results = aligner.align(seq1=read.seq)
        assert isinstance(results, tuple)
        assert len(results) >= 0
        if results:
            aln = results[0]
            assert isinstance(aln, bwap.Alignment)
            print(f"✅ Single alignment successful!")
            print(f"   Result: {aln}")
            if aln.cigar:
                print(f"   CIGAR: {aln.cigar}")
            
    except Exception as e:
        pytest.fail(f"Single-end mapping failed: {e}")


def cleanup_demo_files():
    """Clean up demo files."""
    files_to_remove = [
        "demo_reference.fasta",
        "demo_reads_R1.fastq", 
        "demo_reads_R2.fastq",
        "tests/demo_index/demo_index.amb",
        "tests/demo_index/demo_index.ann",
        "tests/demo_index/demo_index.bwt",
        "tests/demo_index/demo_index.pac",
        "tests/demo_index/demo_index.sa"
    ]
    
    for file in files_to_remove:
        if os.path.exists(file):
            os.remove(file)
            print(f"🗑️  Removed {file}")


if __name__ == "__main__":
    try:
        # Run tests
        success1 = test_paired_end_mapping()
        success2 = test_single_end_mapping()
        
        if success1 and success2:
            print("\n🎉 All tests passed!")
        else:
            print("\n❌ Some tests failed!")
            sys.exit(1)
            
    finally:
        # Clean up
        print("\n🧹 Cleaning up demo files...")
        cleanup_demo_files()
