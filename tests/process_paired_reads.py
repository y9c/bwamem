#!/usr/bin/env python3
"""
Script to process paired-end FASTQ files and map them to a reference genome.

This script:
1. Takes two compressed FASTQ files as input (R1 and R2)
2. Indexes a reference genome
3. Maps read pairs one by one using bwamem
4. Outputs mapping results

Usage:
    python process_paired_reads.py
"""

import os
import sys

# Add the parent directory to the path so we can import bwamem
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from bwamem import BwaAligner, BwaIndexer, read_paired_fastq


def main():
    # Input file paths
    r1_file = os.path.expanduser("~/Desktop/test1.fq.gz")
    r2_file = os.path.expanduser("~/Desktop/test2.fq.gz")
    reference_file = os.path.expanduser("/home/yec/Desktop/genes_rRNA.fa")

    # Check if input files exist
    if not os.path.exists(r1_file):
        print(f"Error: R1 file not found: {r1_file}")
        return 1

    if not os.path.exists(r2_file):
        print(f"Error: R2 file not found: {r2_file}")
        return 1

    if not os.path.exists(reference_file):
        print(f"Error: Reference file not found: {reference_file}")
        return 1

    print("🧬 Paired-End Read Processing with BWA")
    print("=" * 50)
    print(f"R1 file: {r1_file}")
    print(f"R2 file: {r2_file}")
    print(f"Reference: {reference_file}")
    print()

    # Step 1: Build BWA index
    print("🔨 Building BWA index...")
    from pathlib import Path
    test_dir = Path(__file__).parent.parent
    indexer = BwaIndexer()
    index_path = str(test_dir / "test_data" / "rRNA" / "rRNA_index")

    try:
        indexer.build_index(reference_file, index_path)
        print(f"✅ Index built successfully at: {index_path}")
    except Exception as e:
        print(f"❌ Failed to build index: {e}")
        return 1

    # Step 2: Create aligner
    print("🎯 Creating BWA aligner...")
    try:
        aligner = BwaAligner(index_path)
        print("✅ Aligner created successfully")
    except Exception as e:
        print(f"❌ Failed to create aligner: {e}")
        return 1

    # Step 3: Process read pairs
    print("\n📖 Processing read pairs...")
    print("-" * 30)

    total_pairs = 0
    mapped_pairs = 0

    try:
        # Read paired-end FASTQ files
        for read_pair in read_paired_fastq(r1_file, r2_file):
            total_pairs += 1

            # Extract sequences
            seq1 = read_pair[0].sequence
            seq2 = read_pair[1].sequence

            # Align the pair
            try:
                alignments = aligner.align(seq1, seq2)

                if alignments:
                    print(alignments)
                    mapped_pairs += 1
                    print(f"Pair {total_pairs}: {len(alignments)} alignment(s)")

                    # Print details for each alignment
                    for i, aln in enumerate(alignments):
                        if hasattr(aln, "read1") and hasattr(aln, "read2"):
                            # Paired alignment
                            print(f"  Alignment {i + 1}:")
                            print(
                                f"    Read1: {aln.read1.rname}:{aln.read1.pos} {aln.read1.orient} "
                                f"mapq={aln.read1.mapq} cigar={aln.read1.cigar}"
                            )
                            print(
                                f"    Read2: {aln.read2.rname}:{aln.read2.pos} {aln.read2.orient} "
                                f"mapq={aln.read2.mapq} cigar={aln.read2.cigar}"
                            )
                            print(f"    Proper pair: {aln.is_proper_pair}")
                            if aln.insert_size is not None:
                                print(f"    Insert size: {aln.insert_size}")
                        else:
                            # Single alignment
                            print(
                                f"  Alignment {i + 1}: {aln.rname}:{aln.pos} {aln.orient} "
                                f"mapq={aln.mapq} cigar={aln.cigar}"
                            )
                else:
                    print(f"Pair {total_pairs}: No alignments found")

            except Exception as e:
                print(f"Pair {total_pairs}: Error during alignment: {e}")

            # Print progress every 100 pairs
            if total_pairs % 100 == 0:
                print(f"Processed {total_pairs} pairs...")

    except Exception as e:
        print(f"❌ Error processing reads: {e}")
        return 1

    # Step 4: Summary
    print("\n📊 Summary")
    print("=" * 20)
    print(f"Total read pairs processed: {total_pairs}")
    print(f"Pairs with alignments: {mapped_pairs}")
    if total_pairs > 0:
        mapping_rate = (mapped_pairs / total_pairs) * 100
        print(f"Mapping rate: {mapping_rate:.1f}%")

    print("\n✅ Processing completed successfully!")

    # Clean up index files (optional)
    cleanup = input("\n🗑️  Remove index files? (y/N): ").strip().lower()
    if cleanup in ["y", "yes"]:
        try:
            for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
                index_file = f"{index_path}{ext}"
                if os.path.exists(index_file):
                    os.remove(index_file)
            print("✅ Index files cleaned up")
        except Exception as e:
            print(f"⚠️  Warning: Could not clean up some index files: {e}")

    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
