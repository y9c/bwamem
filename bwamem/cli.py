#!/usr/bin/env python3
"""bwamem — BWA-MEM Python bindings CLI.

Usage:
    bwamem map -i ref -k 18 -n 0.05 reads.fq
    bwamem map -i cont -i gene -i tx -k 18,14,10 -n 0.05,0.15,0.30 reads.fq
    bwamem map -i ref -1 r1.fq -2 r2.fq
    bwamem index ref1.fa ref2.fa ref3.fa
"""

import argparse
import os
import sys
import time

from bwamem import HierarchicalAligner, BwaIndexer, FastxReader


def _parse_params(args):
    refs = args.index
    ks = [int(x) for x in args.seed_len.split(",")] if args.seed_len else [19]
    ns = [float(x) for x in args.nm_ratio.split(",")] if args.nm_ratio else [1.0]
    ss = [int(x) for x in args.min_score.split(",")] if args.min_score else [30]
    if len(ks) == 1: ks = ks * len(refs)
    if len(ns) == 1: ns = ns * len(refs)
    if len(ss) == 1: ss = ss * len(refs)
    if len(ks) != len(refs) or len(ns) != len(refs) or len(ss) != len(refs):
        print("error: -k/-n/-T must have 1 or N values (N=number of -i refs)", file=sys.stderr)
        sys.exit(1)
    layers = []
    for i, ref in enumerate(refs):
        layers.append({
            "index_prefix": ref, "min_seed_len": ks[i],
            "max_nm_ratio": ns[i], "min_score": ss[i],
        })
    return layers


def _sam_header(layers, rg_id=None, rg_sm=None):
    lines = ["@HD\tVN:1.6\tSO:queryname"]
    for i, layer in enumerate(layers):
        a = layer["aligner"]
        for j, name in enumerate(a._rid_to_name):
            length = 0
            try:
                seq_str = a.seq_by_rid(j, 0, 1)
                if seq_str is not None:
                    length = len(seq_str)
            except Exception:
                pass
            lines.append(f"@SQ\tSN:{name}\tLN:{length}")
    if rg_id or rg_sm:
        rg_line = "@RG"
        if rg_id:
            rg_line += f"\tID:{rg_id}"
        if rg_sm:
            rg_line += f"\tSM:{rg_sm}"
        lines.append(rg_line)
    return "\n".join(lines)


def _se_layer_name(i, layers):
    return os.path.basename(layers[i]["index_prefix"])


def _pe_layer_name(i, layers):
    prefix = layers[i]["index_prefix"]
    return os.path.basename(prefix)


def cmd_map(args):
    layers_cfg = _parse_params(args)
    aligner = HierarchicalAligner(layers_cfg)
    pe_mode = bool(args.r1)

    out_fh = open(args.output, "w") if args.output else sys.stdout
    try:
        header = _sam_header(aligner._layers, rg_id=args.rg_id, rg_sm=args.rg_sm)
        out_fh.write(header + "\n")
        n_layers = aligner.n_layers
        layer_counts = [0] * n_layers
        n_total = 0
        t0 = time.time()

        if pe_mode:
            reader = FastxReader(args.r1, args.r2)
            for pair in reader:
                r1, r2 = pair
                results, layer = aligner.align_pe(r1.sequence, r2.sequence, min_mapq=0)
                n_total += 1
                packed = _pack_pe_results(results, layer, r1, r2)
                if layer >= 0:
                    layer_counts[layer] += 1
                out_fh.write(packed + "\n")
        else:
            reader = FastxReader(args.reads)
            for read in reader:
                hits, layer = aligner.align(read.sequence, min_mapq=0)
                n_total += 1
                if layer >= 0:
                    layer_counts[layer] += 1
                qual = read.quality if read.quality else "*"
                flag = 0 if layer >= 0 else 4
                rname = hits[0][0] if hits else "*"
                pos = int(hits[0][1]) + 1 if hits else 0
                mapq = hits[0][6] if hits else 0
                cigar = hits[0][7] if hits else "*"
                hi_tag = f"HI:Z:{layer}" if layer >= 0 else "HI:Z:-1"
                out_fh.write(
                    f"{read.name}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}"
                    f"\t*\t0\t0\t{read.sequence}\t{qual}\t{hi_tag}\n"
                )
    finally:
        if args.output:
            out_fh.close()

    elapsed = time.time() - t0
    msg = f"[bwamem] {n_total} reads, {n_total * len(layer_counts)} lines in {elapsed:.1f}s"
    print(msg, file=sys.stderr)
    if args.report:
        _write_report(args.report, n_total, layer_counts, layers_cfg)


def _pack_pe_results(results, layer, r1, r2):
    """Pack PE results into two SAM lines (R1 then R2)."""
    qual1 = r1.quality if r1.quality else "*"
    qual2 = r2.quality if r2.quality else "*"
    hi = f"HI:Z:{layer}" if layer >= 0 else "HI:Z:-1"
    if layer >= 0 and results:
        h1, h2, is_paired, isize = results[0]
        if h1 is None or h2 is None or not is_paired:
            pairs = [(h1, r1, qual1, False), (h2, r2, qual2, True)]
        else:
            pairs = [(h1, r1, qual1, False), (h2, r2, qual2, True)]
    else:
        pairs = [(None, r1, qual1, False), (None, r2, qual2, True)]

    lines = []
    for h, read, qual, is_r2 in pairs:
        if h is None:
            flag = 77 if is_r2 else 69
            lines.append(f"{read.name}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t{read.sequence}\t{qual}\t{hi}")
        else:
            flag = 83 if is_r2 else 99
            mpos = int(h[1]) + 1
            isize_val = abs(results[0][3]) if results and results[0][2] else 0
            lines.append(
                f"{read.name}\t{flag}\t{h[0]}\t{mpos}\t{h[6]}\t{h[7]}"
                f"\t=\t{mpos}\t{isize_val}\t{read.sequence}\t{qual}\t{hi}"
            )
    return "\n".join(lines)


def _write_report(path, n_total, layer_counts, layers_cfg):
    with open(path, "w") as f:
        f.write(f"Total reads (pairs): {n_total}\n")
        for i, c in enumerate(layer_counts):
            name = os.path.basename(layers_cfg[i]["index_prefix"])
            pct = c / n_total * 100 if n_total else 0
            f.write(f"Layer {i} ({name}): {c} ({pct:.1f}%)\n")
        unmapped = n_total - sum(layer_counts)
        f.write(f"Unmapped: {unmapped} ({unmapped/n_total*100:.1f}%)\n" if n_total else "Unmapped: 0\n")


def cmd_index(args):
    indexer = BwaIndexer()
    for fa in args.fasta:
        prefix = os.path.splitext(fa)[0]
        print(f"[bwamem] indexing {fa} -> {prefix}.*", file=sys.stderr)
        indexer.build_index(fa, prefix)
    print("[bwamem] indexing done", file=sys.stderr)


def cmd_run_layer(args):
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
                       help="min seed length, comma-separated per index")
    p_map.add_argument("-n", "--nm-ratio", default="1.0",
                       help="max NM/len ratio, comma-separated per index")
    p_map.add_argument("-T", "--min-score", default="30",
                       help="min alignment score, comma-sep per index")
    p_map.add_argument("-o", "--output", default=None,
                       help="output SAM file (default: stdout)")
    p_map.add_argument("--rg-id", default=None, help="read group ID")
    p_map.add_argument("--rg-sm", default=None, help="read group sample name")
    p_map.add_argument("--report", default=None, help="write per-layer stats to file")
    p_map.add_argument("-1", "--r1", default=None, help="R1 FASTQ (PE mode)")
    p_map.add_argument("-2", "--r2", default=None, help="R2 FASTQ (PE mode)")
    p_map.add_argument("reads", nargs="?", default=None, help="FASTQ file (SE mode)")
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
    if args.cmd == "map" and not args.r1 and not args.r2 and not args.reads:
        print("error: specify reads (positional for SE, -1/-2 for PE)", file=sys.stderr)
        return 1
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
