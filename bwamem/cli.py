#!/usr/bin/env python3
"""bwamem — BWA-MEM Python bindings CLI.

Usage:
    bwamem map -i ref -k 18 -n 0.05 reads.fq
    bwamem map -i cont -i gene -i tx -k 18,14,10 -n 0.05,0.15,0.30 reads.fq
    bwamem index ref1.fa ref2.fa ref3.fa
"""

import argparse
import os
import sys
import time

from bwamem import HierarchicalAligner, BwaIndexer, FastxReader



def cmd_map(args):
    refs = args.index
    ks = [int(x) for x in args.seed_len.split(",")] if args.seed_len else [19]
    ns = [float(x) for x in args.nm_ratio.split(",")] if args.nm_ratio else [1.0]
    ss = [int(x) for x in args.min_score.split(",")] if args.min_score else [30]

    if len(ks) == 1: ks = ks * len(refs)
    if len(ns) == 1: ns = ns * len(refs)
    if len(ss) == 1: ss = ss * len(refs)

    if len(ks) != len(refs) or len(ns) != len(refs) or len(ss) != len(refs):
        print("error: --seed-len/--nm-ratio/--min-score must have 1 or N values (N=number of -i refs)", file=sys.stderr)
        return 1

    layers = []
    for i, ref in enumerate(refs):
        layers.append({
            "index_prefix": ref, "min_seed_len": ks[i],
            "max_nm_ratio": ns[i], "min_score": ss[i],
        })
    aligner = HierarchicalAligner(layers)
    reader = FastxReader(args.reads)
    t0 = time.time()
    n_total, n_mapped = 0, 0
    layer_counts = [0] * len(layers)

    for read in reader:
        hits, layer = aligner.align(read.sequence, min_mapq=0)
        n_total += 1
        if layer >= 0:
            n_mapped += 1
            layer_counts[layer] += 1

        # emit SAM line
        flag = 0 if layer >= 0 else 4
        rname = "*"
        pos = 0
        mapq = 0
        cigar = "*"
        if hits:
            rname = hits[0][0]       # contig name
            pos = int(hits[0][1]) + 1  # 1-based
            mapq = hits[0][6]
            cigar = hits[0][7]

        hi_tag = f"HI:Z:{layer}" if layer >= 0 else "HI:Z:-1"
        qual = read.quality if read.quality else "*"
        print(f"{read.name}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{read.sequence}\t{qual}\t{hi_tag}")

    elapsed = time.time() - t0
    print(f"[bwamem] {n_total} reads, {n_mapped} mapped "
          f"({','.join(f'L{i}={layer_counts[i]}' for i in range(len(layers)))}), "
          f"{elapsed:.1f}s", file=sys.stderr)


def cmd_index(args):
    indexer = BwaIndexer()
    for fa in args.fasta:
        prefix = os.path.splitext(fa)[0]
        print(f"[bwamem] indexing {fa} → {prefix}.*", file=sys.stderr)
        indexer.build_index(fa, prefix)
    print("[bwamem] indexing done", file=sys.stderr)


def cmd_run_layer(args):
    """Run a single layer (non-hierarchical) — useful for debugging."""
    from bwamem import BwaAligner

    aligner = BwaAligner(
        index_prefix=args.index,
        min_seed_len=args.min_seed_len,
        min_score=args.min_score,
    )
    reader = FastxReader(args.reads)
    for read in reader:
        hits = aligner.align(read.sequence, min_mapq=0)
        flag = 0 if hits else 4
        rname = hits[0][0] if hits else "*"
        pos = int(hits[0][1]) + 1 if hits else 0
        mapq = hits[0][6] if hits else 0
        cigar = hits[0][7] if hits else "*"
        qual = read.quality if read.quality else "*"
        print(f"{read.name}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{read.sequence}\t{qual}")


def main():
    ap = argparse.ArgumentParser(description="bwamem: BWA-MEM Python bindings CLI")
    sub = ap.add_subparsers(dest="cmd")

    p_map = sub.add_parser("map", help="BWA-MEM mapping (single or hierarchical)")
    p_map.add_argument("-i", "--index", action="append", required=True,
                       help="BWA index prefix (repeatable)")
    p_map.add_argument("-k", "--seed-len", default="19",
                       help="min seed length, comma-separated per index (default: 19)")
    p_map.add_argument("-n", "--nm-ratio", default="1.0",
                       help="max NM/len ratio, comma-separated per index (default: 1.0)")
    p_map.add_argument("-T", "--min-score", default="30",
                       help="min alignment score, comma-sep per index (default: 30)")
    p_map.add_argument("reads", help="FASTQ file (or - for stdin)")
    p_map.set_defaults(func=cmd_map)

    p_idx = sub.add_parser("index", help="build BWA indices")
    p_idx.add_argument("fasta", nargs="+", help="FASTA file(s)")
    p_idx.set_defaults(func=cmd_index)

    p_layer = sub.add_parser("layer", help="single-layer mapping (debug)")
    p_layer.add_argument("-i", "--index", required=True)
    p_layer.add_argument("-k", "--min-seed-len", type=int, default=19)
    p_layer.add_argument("-T", "--min-score", type=int, default=30)
    p_layer.add_argument("reads")
    p_layer.set_defaults(func=cmd_run_layer)

    args = ap.parse_args()
    if args.cmd is None:
        ap.print_help()
        return 1
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
