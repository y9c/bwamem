#!/usr/bin/env python3
"""
Simple script to demonstrate paired-end read mapping with bwamem.

This script shows how to:
1. Index a reference genome
2. Process paired-end FASTQ files
3. Map read pairs and display results

Usage:
    python simple_paired_mapping.py
"""

import os
import sys
import gzip
from pathlib import Path

# Add the parent directory to the path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from bwamem import BwaAligner, BwaIndexer, read_paired_fastq
except ImportError as e:
    print(f"Error importing bwamem: {e}")
    print("Make sure you're running this from the bwamem project directory")
    sys.exit(1)


def process_paired_reads():
    """Process paired-end reads and map them to reference."""
    
    # File paths
    r1_file = os.path.expanduser("~/Desktop/test1.fq.gz")
    r2_file = os.path.expanduser("~/Desktop/test2.fq.gz")
    reference_file = os.path.expanduser("/home/yec/Desktop/genes_rRNA.fa")
    index_path = "rRNA_index"
    
    print("🧬 BWA Paired-End Read Mapping Demo")
    print("=" * 40)
    
    # Check files exist
    for file_path, name in [(r1_file, "R1"), (r2_file, "R2"), (reference_file, "Reference")]:
        if not os.path.exists(file_path):
            print(f"❌ {name} file not found: {file_path}")
            return False
        print(f"✅ {name}: {file_path}")
    
    print()
    
    # Build index
    print("🔨 Building BWA index...")
    try:
        indexer = BwaIndexer()
        indexer.build_index(reference_file, index_path)
        print(f"✅ Index built: {index_path}")
    except Exception as e:
        print(f"❌ Index building failed: {e}")
        return False
    
    # Create aligner
    print("🎯 Creating aligner...")
    try:
        aligner = BwaAligner(index_path)
        print("✅ Aligner ready")
    except Exception as e:
        print(f"❌ Aligner creation failed: {e}")
        return False
    
    # Process reads
    print("\n📖 Processing read pairs...")
    print("-" * 30)
    
    pair_count = 0
    mapped_count = 0
    
    try:
        for read_pair in read_paired_fastq(r1_file, r2_file):
            pair_count += 1
            
            # Get sequences
            seq1 = read_pair.read1.sequence
            seq2 = read_pair.read2.sequence
            
            print(f"Pair {pair_count}: {read_pair.read1.name}")
            
            # Align
            try:
                alignments = aligner.align(seq1, seq2)
                
                if alignments:
                    mapped_count += 1
                    print(f"  ✅ {len(alignments)} alignment(s) found")
                    
                    for i, aln in enumerate(alignments):
                        if hasattr(aln, 'read1') and hasattr(aln, 'read2'):
                            # Paired alignment
                            print(f"    Alignment {i+1}:")
                            print(f"      R1: {aln.read1.chromosome}:{aln.read1.position} "
                                  f"{aln.read1.strand} (mapq={aln.read1.mapq})")
                            print(f"      R2: {aln.read2.chromosome}:{aln.read2.position} "
                                  f"{aln.read2.strand} (mapq={aln.read2.mapq})")
                            print(f"      Proper pair: {aln.is_proper_pair}")
                        else:
                            # Single alignment
                            print(f"    Alignment {i+1}: {aln.chromosome}:{aln.position} "
                                  f"{aln.strand} (mapq={aln.mapq})")
                else:
                    print("  ❌ No alignments found")
                    
            except Exception as e:
                print(f"  ⚠️  Alignment error: {e}")
            
            # Limit to first 10 pairs for demo
            if pair_count >= 10:
                print(f"\n... (showing first {pair_count} pairs only)")
                break
                
    except Exception as e:
        print(f"❌ Read processing error: {e}")
        return False
    
    # Summary
    print(f"\n📊 Results:")
    print(f"  Processed: {pair_count} pairs")
    print(f"  Mapped: {mapped_count} pairs")
    if pair_count > 0:
        rate = (mapped_count / pair_count) * 100
        print(f"  Mapping rate: {rate:.1f}%")
    
    print("\n✅ Demo completed!")
    return True


if __name__ == "__main__":
    success = process_paired_reads()
    sys.exit(0 if success else 1)