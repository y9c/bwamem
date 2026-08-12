#!/usr/bin/env python3
"""bwamem-hier — hierarchical multi-reference BWA-MEM mapper (CLI).

Usage:
    bwamem-hier map -r cont,k=18,n=0.05 -r gene,k=14,n=0.15 \
                    -r tx,k=10,n=0.30 reads.fq > out.sam

    bwamem-hier index ref1.fa ref2.fa ref3.fa
"""

import argparse
import os
import sys
import time

from bwamem import HierarchicalAligner, BwaIndexer, FastxReader


def parse_layer(spec: str) -> dict:
    cfg = {}
    for part in spec.split(","):
        k, _, v = part.partition("=")
        k = k.strip()
        v = v.strip()
        if k in ("path", "p", "index"):
            cfg["index_prefix"] = v
        elif k in ("k", "min_seed_len"):
            cfg["min_seed_len"] = int(v)
        elif k in ("n", "nm", "max_nm_ratio"):
            cfg["max_nm_ratio"] = float(v)
        elif k in ("s", "min_score"):
            cfg["min_score"] = int(v)
        elif k in ("c", "max_occ"):
            cfg["max_occ"] = int(v)
    return cfg


def cmd_map(args):
    layers = [parse_layer(s) for s in args.ref]
    if not layers:
        print("error: at least one -r/--ref required", file=sys.stderr)
        return 1

    aligner = HierarchicalAligner(layers)
    reader = FastxReader(args.reads)
    t0 = time.time()
    n_total, n_mapped = 0, 0
    layer_counts = [0] * len(layers)

    for name, seq, qual in reader:
        hits, layer = aligner.align(seq, min_mapq=0)
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
        print(f"{name}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\t{hi_tag}")

    elapsed = time.time() - t0
    print(f"[bwamem-hier] {n_total} reads, {n_mapped} mapped "
          f"({','.join(f'L{i}={layer_counts[i]}' for i in range(len(layers)))}), "
          f"{elapsed:.1f}s", file=sys.stderr)


def cmd_index(args):
    indexer = BwaIndexer()
    for fa in args.fasta:
        prefix = os.path.splitext(fa)[0]
        print(f"[bwamem-hier] indexing {fa} → {prefix}.*", file=sys.stderr)
        indexer.build_index(fa, prefix)
    print("[bwamem-hier] indexing done", file=sys.stderr)


def cmd_run_layer(args):
    """Run a single layer (non-hierarchical) — useful for debugging."""
    from bwamem import BwaAligner

    aligner = BwaAligner(
        index_prefix=args.index,
        min_seed_len=args.min_seed_len,
        min_score=args.min_score,
    )
    reader = FastxReader(args.reads)
    for name, seq, qual in reader:
        hits = aligner.align(seq, min_mapq=0)
        flag = 0 if hits else 4
        rname = hits[0][0] if hits else "*"
        pos = int(hits[0][1]) + 1 if hits else 0
        mapq = hits[0][6] if hits else 0
        cigar = hits[0][7] if hits else "*"
        print(f"{name}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}")


def main():
    ap = argparse.ArgumentParser(description="bwamem-hier: hierarchical BWA-MEM mapper")
    sub = ap.add_subparsers(dest="cmd")

    p_map = sub.add_parser("map", help="hierarchical mapping")
    p_map.add_argument("-r", "--ref", action="append", required=True,
                       help="layer spec: path,k=SEED,n=NM_RATIO (repeatable)")
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
