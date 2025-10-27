#!/usr/bin/env python3
"""
Example usage of the bwap package with both single-end and paired-end alignment support.
"""

from bwap import BwaAligner, BwaIndexer, Alignment, PairedAlignment
import os

def main():
    """Demonstrate basic usage of bwap with SE and PE alignment."""
    print("Bwap - Python bindings for BWA aligner")
    print("=====================================")
    print("Supports both Single-End (SE) and Paired-End (PE) alignment")
    
    # Demonstrate BwaIndexer functionality
    print("\n1. BWA Index Building:")
    print("---------------------")
    
    indexer = BwaIndexer()
    print(f"Available algorithms: {['auto', 'rb2', 'bwtsw', 'is']}")
    print(f"Default algorithm: {indexer.algorithm}")
    print(f"Default block size: {indexer.block_size}")
    
    # Example of how to build an index (would need actual FASTA file)
    print("\nTo build an index from a FASTA file:")
    print("indexer = BwaIndexer(algorithm='auto')")
    print("index_path = indexer.build_index('reference.fa')")
    print("print(f'Index built at: {index_path}')")
    
    # Demonstrate BwaAligner functionality
    print("\n2. BWA Alignment:")
    print("----------------")
    
    try:
        # This will fail without a valid index, but shows the API
        aligner = BwaAligner("path/to/your/bwa/index")
        print("BwaAligner initialized successfully!")
        
        # Example sequences
        seq1 = "ACGATCGCGATCGA"
        seq2 = "TTCGATCGATCGAT"
        print(f"Example sequences:")
        print(f"  Read 1: {seq1}")
        print(f"  Read 2: {seq2}")
        
    except ValueError as e:
        print(f"Expected error (no valid index): {e}")
        print("To use this package, you need to:")
        print("1. Build a BWA index: indexer.build_index('reference.fa')")
        print("2. Use the index path when creating BwaAligner")
    
    print("\n3. Single-End vs Paired-End Alignment:")
    print("--------------------------------------")
    print("The new align() method supports both single-end and paired-end reads:")
    print("")
    print("# Single-end alignment:")
    print("alignments = aligner.align('ACGATCGCGATCGA')")
    print("# Returns tuple of Alignment objects")
    print("for aln in alignments:")
    print("    print(f'  {aln.rname}:{aln.pos} {aln.orient} {aln.cigar} (mapq={aln.mapq}, score={aln.score})')")
    print("")
    print("# Paired-end alignment:")
    print("paired_alignments = aligner.align('ACGATCGCGATCGA', 'TTCGATCGATCGAT')")
    print("# Returns tuple of PairedAlignment objects")
    print("for pe_aln in paired_alignments:")
    print("    print(f'  Read1: {pe_aln.read1.rname}:{pe_aln.read1.pos} {pe_aln.read1.orient}')")
    print("    print(f'  Read2: {pe_aln.read2.rname}:{pe_aln.read2.pos} {pe_aln.read2.orient}')")
    print("    print(f'  Proper pair: {pe_aln.is_proper_pair}, Insert size: {pe_aln.insert_size}')")
    print("")
    print("# Paired-end with custom insert size:")
    print("paired_alignments = aligner.align('ACGATCGCGATCGA', 'TTCGATCGATCGAT',")
    print("                                 insert_size=500, insert_std=50)")
    print("")
    print("# Data structures:")
    print("Alignment: (rname, orient, pos, mapq, cigar, NM, score, is_primary)")
    print("PairedAlignment: (read1, read2, is_proper_pair, insert_size)")
    
    print("\n4. Complete Workflow Example:")
    print("----------------------------")
    print("# Step 1: Build index from FASTA")
    print("indexer = BwaIndexer(algorithm='auto')")
    print("index_path = indexer.build_index('reference.fa')")
    print("")
    print("# Step 2: Create aligner with the index")
    print("aligner = BwaAligner(index_path)")
    print("")
    print("# Step 3a: Single-end alignment")
    print("se_alignments = aligner.align('ACGATCGCGATCGA')")
    print("for aln in se_alignments:")
    print("    print(f'  {aln.rname}:{aln.pos} {aln.orient} {aln.cigar} (mapq={aln.mapq})')")
    print("")
    print("# Step 3b: Paired-end alignment")
    print("pe_alignments = aligner.align('ACGATCGCGATCGA', 'TTCGATCGATCGAT')")
    print("for pe_aln in pe_alignments:")
    print("    print(f'  Read1: {pe_aln.read1.rname}:{pe_aln.read1.pos}')")
    print("    print(f'  Read2: {pe_aln.read2.rname}:{pe_aln.read2.pos}')")
    print("    print(f'  Proper pair: {pe_aln.is_proper_pair}, Insert size: {pe_aln.insert_size}')")
    
    print("\n5. Advanced Paired-End Features:")
    print("--------------------------------")
    print("# Custom insert size distribution")
    print("pe_alignments = aligner.align('read1.fastq', 'read2.fastq',")
    print("                             insert_size=300, insert_std=30)")
    print("")
    print("# Filter for proper pairs only")
    print("proper_pairs = [pe for pe in pe_alignments if pe.is_proper_pair]")
    print("")
    print("# Access individual read alignments")
    print("for pe_aln in pe_alignments:")
    print("    read1_aln = pe_aln.read1")
    print("    read2_aln = pe_aln.read2")
    print("    if read1_aln.is_primary and read2_aln.is_primary:")
    print("        print(f'Primary alignment: {read1_aln.rname}:{read1_aln.pos}-{read2_aln.pos}')")
    
    print("\nPackage information:")
    print(f"Version: {BwaAligner.__module__}")

if __name__ == "__main__":
    main()
