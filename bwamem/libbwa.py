import argparse
import importlib
import importlib.util
import os
import re
import sys
import threading
from collections import deque, namedtuple
from contextlib import contextmanager
from typing import Optional, Tuple

from cffi import FFI

# CFFI configuration
ffi = FFI()


def get_shared_lib(name):
    """Locate and open the shared library."""
    lib_file = None
    current_dir = os.path.dirname(os.path.abspath(__file__))
    for f in os.listdir(current_dir):
        if f.startswith(name) and (
            f.endswith(".so") or f.endswith(".dll") or f.endswith(".dylib")
        ):
            lib_file = os.path.join(current_dir, f)
            break
    if lib_file is None:
        try:
            spec = importlib.util.find_spec(name)
            if spec is None or spec.origin is None:
                spec = importlib.util.find_spec("bwamem." + name)
            if spec is None or spec.origin is None:
                raise ImportError(f'Cannot locate C library "{name}".')
            lib_file = spec.origin
        except Exception:
            raise ImportError('Cannot locate C library "{}".'.format(name))
    lib_file = os.path.abspath(lib_file)
    library = ffi.dlopen(lib_file)
    return library


libbwa = get_shared_lib("bwalib")

ffi.cdef("""
  typedef struct {
    int64_t rb, re;
    int qb, qe;
    int rid;
    int score;
    int truesc;
    int sub;
    int alt_sc;
    int csub;
    int sub_n;
    int w;
    int seedcov;
    int secondary;
    int secondary_all;
    int seedlen0;
    int n_comp:30, is_alt:2;
    float frac_rep;
    uint64_t hash;
  } mem_alnreg_t;

  typedef struct { size_t n, m; mem_alnreg_t *a; } mem_alnreg_v;

  typedef struct {
    int64_t pos;
    int rid;
    int flag;
    uint32_t is_rev:1, is_alt:1, mapq:8, NM:22;
    int n_cigar;
    uint32_t *cigar;
    int score;
    int sub;
    int alt_sc;
  } mem_aln_t;

  typedef struct {
    int a, b, q, e, q2, e2;
    int pen_unpaired;
    int pen_clip5, pen_clip3;
    int w, zdrop, max_chain_extend, min_seed_len;
    int T, flag, min_chain_weight;
    int max_occ, max_chain_gap;
    int n_threads, chunk_size;
    float mask_level, chain_drop_ratio, mask_level_red;
    float Gdrop_frac;
    float split_factor;
    int split_width;
    int max_matesw;
    int mapQ_coef_len, mapQ_coef_fac;
    int max_ins;
    int rev_len2;
    int dump_prefix;
  } mem_opt_t;

  typedef struct {
    int64_t offset;
    int32_t len;
    int32_t n_amb, gi;
    char *name, *anno;
  } bntann1_t;

  typedef struct {
    int64_t l_pac;
    int32_t n_seqs;
    uint32_t seed;
    bntann1_t *anns;
    int32_t n_holes;
    void *ambs;
    FILE *fp_pac;
  } bntseq_t;

  typedef struct {
    void *bwt;
    bntseq_t *bns;
    uint8_t *pac;
  } bwaidx_t;

  typedef struct {
    int low, high, failed;
    double avg, std;
  } mem_pestat_t;

  typedef struct {
    int l_seq, id;
    char *name, *comment, *seq, *qual, *sam;
  } bseq1_t;

  typedef struct {
    size_t l, m;
    char *s;
  } kstring_t;

  typedef struct {
    kstring_t name, comment, seq, qual;
    int last_char;
    int f;
  } kseq_t;

  mem_opt_t *mem_opt_init();
  void bwa_fill_scmat(int a, int b, int8_t mat[25]);
  bwaidx_t *bwa_idx_load(const char *hint, int which);
  void bwa_idx_destroy(bwaidx_t *idx);
  mem_alnreg_v mem_align1(const mem_opt_t *opt, const void *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq);
  mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const char *query_, const mem_alnreg_t *ar);
  mem_opt_t *get_opts(int argc, char *argv[], const bwaidx_t *idx);
  mem_aln_t *mem_reg2aln_ptr(const mem_opt_t* opt, const bntseq_t* bns, const uint8_t* pac, int l_seq, const char* seq, const mem_alnreg_t* ar);
  void mem_pestat(const mem_opt_t *opt, int64_t l_pac, int n, const mem_alnreg_v *regs, mem_pestat_t pes[4]);
  void mem_matesw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], const mem_alnreg_t *a, int l_ms, const uint8_t *ms, mem_alnreg_v *regs);
  int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id);
  uint8_t *encode_seq(const char *seq, int len);
  char *build_cigar_string(const uint32_t *cigar, int n_cigar);
  void free(void *ptr);

  void *err_xzopen_core(const char *prog, const char *fn, const char *mode);
  void gzclose(void *fp);
  kseq_t *kseq_init(void *fp);
  void kseq_destroy(kseq_t *ks);
  int kseq_read(kseq_t *seq);
  bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1, void *ks2);

  typedef struct { uint32_t len; uint32_t op; } cigar_pair_t;
  int64_t compute_r_en(int64_t r_st, const cigar_pair_t* cigar, int n_cigar);
  int compute_blen(const cigar_pair_t* cigar, int n_cigar);
  int compute_mlen(const cigar_pair_t* cigar, int n_cigar);
  char *bwa_fetch_seq(const bwaidx_t *idx, int rid, int64_t start, int64_t end);
  
  typedef struct {
    char* cigar;
    long long pos;
    long long r_en;
    int rid;
    int strand;
    int q_st;
    int q_en;
    int mapq;
    int NM;
    int score;
  } raw_hit_t;

  typedef struct {
    size_t n;
    raw_hit_t* hits;
  } raw_hit_v;

  typedef struct {
    mem_alnreg_v regs[2];
    mem_pestat_t pes[4];
  } pe_regs_t;

  raw_hit_v* bwa_mem_reg2aln_all(const mem_opt_t* opt, const bntseq_t* bns, const uint8_t* pac, int l_seq, const char* seq, const mem_alnreg_v* regs, int min_mapq, int min_blen, int min_mlen);
  void free_raw_hit_v(raw_hit_v* v);

  pe_regs_t* bwa_align_pe(const mem_opt_t* opt, const bwaidx_t* idx, const char* seq1, int l1, const char* seq2, int l2, double avg, double std, int imin, int imax);
  void free_pe_regs(pe_regs_t* p);

  char* visualize_alignment_c(const char* ctg, int64_t r_st, int64_t r_en, int strand, int score, int mapq, const char* ref_seq, int ref_len, const char* query_seq, int query_len, const cigar_pair_t* cigar, int n_cigar, int q_st, int q_en, int lw);

  int bwa_idx_build(const char *fa, const char *prefix, int algo_type, int block_size);
  extern int bwa_verbose;
""")


@contextmanager
def suppress_stderr():
    old_stderr_fd = os.dup(2)
    try:
        with open(os.devnull, "w") as devnull:
            os.dup2(devnull.fileno(), 2)
            yield
    finally:
        os.dup2(old_stderr_fd, 2)
        os.close(old_stderr_fd)


class BwaAligner:
    def __init__(
        self,
        index_prefix,
        min_seed_len=19,
        max_occ=500,
        softclip_supplementary=True,
        mark_secondary=True,
        clip_penalties=(6, 6),
        unpaired_penalty=24,
        min_score=30,
        insert_model=None,
    ):
        try:
            libbwa.bwa_verbose = 0
        except:
            pass
        with suppress_stderr():
            self.index = libbwa.bwa_idx_load(index_prefix.encode(), 7)
        if self.index == ffi.NULL:
            raise ValueError(f"Failed to load BWA index: {index_prefix}")
        argv = [
            "bwamem",
            "-k",
            str(min_seed_len),
            "-c",
            str(max_occ),
            "-L",
            f"{clip_penalties[0]},{clip_penalties[1]}",
            "-U",
            str(unpaired_penalty),
            "-T",
            str(min_score),
        ]
        if softclip_supplementary:
            argv.append("-Y")
        if mark_secondary:
            argv.append("-M")
        argv.append(index_prefix)
        c_strs = [ffi.new("char[]", x.encode()) for x in argv]
        # Create a NULL-terminated array of pointers and keep it alive as an attribute
        self._c_argv = ffi.new("char*[]", len(argv) + 1)
        for i, s in enumerate(c_strs):
            self._c_argv[i] = s
        self._c_argv[len(argv)] = ffi.NULL
        # Store strings too to prevent GC
        self._c_strs = c_strs

        self.opt = libbwa.get_opts(len(argv), self._c_argv, self.index)
        if self.opt == ffi.NULL:
            raise RuntimeError("Failed to initialize BWA options")

        self._insert_model = insert_model
        self._rid_to_name = [
            ffi.string(self.index.bns.anns[i].name).decode()
            for i in range(self.index.bns.n_seqs)
        ]
        self._name_to_rid = {name: i for i, name in enumerate(self._rid_to_name)}

    def __del__(self):
        if hasattr(self, "index") and self.index != ffi.NULL:
            libbwa.bwa_idx_destroy(self.index)
        if hasattr(self, "opt") and self.opt != ffi.NULL:
            libbwa.free(self.opt)

    def seq(self, name, start=0, end=0x7FFFFFFF):
        rid = self._name_to_rid.get(name, -1)
        if rid < 0:
            return None
        return self.seq_by_rid(rid, start, end)

    def seq_by_rid(self, rid, start=0, end=0x7FFFFFFF):
        if rid < 0 or rid >= len(self._rid_to_name):
            return None
        res_ptr = libbwa.bwa_fetch_seq(self.index, rid, start, end)
        if res_ptr == ffi.NULL:
            return None
        res = ffi.string(res_ptr).decode()
        libbwa.free(res_ptr)
        return res

    def align(self, seq1, seq2=None, min_mapq=0, min_blen=0, min_mlen=0):
        if seq2 is None:
            return self.align_raw(seq1, min_mapq, min_blen, min_mlen)
        return self.align_raw_pe(seq1, seq2, min_mapq, min_blen, min_mlen)

    def align_raw(self, seq, min_mapq=0, min_blen=0, min_mlen=0):
        regs = libbwa.mem_align1(
            self.opt,
            self.index.bwt,
            self.index.bns,
            self.index.pac,
            len(seq),
            seq.encode(),
        )
        if regs.n == 0:
            if regs.a != ffi.NULL:
                libbwa.free(regs.a)
            return []
        hits = self._conv_hits_raw(regs, seq, min_mapq, min_blen, min_mlen)
        if regs.a != ffi.NULL:
            libbwa.free(regs.a)
        return hits

    def align_raw_pe(self, seq1, seq2, min_mapq=0, min_blen=0, min_mlen=0):
        model = self._insert_model
        avg = model[0] if model else 0.0
        std = model[1] if model and len(model) > 1 else (avg * 0.1 if avg else 0.0)
        imax = model[2] if model and len(model) > 2 else 0
        imin = model[3] if model and len(model) > 3 else 0
        with suppress_stderr():
            pe = libbwa.bwa_align_pe(
                self.opt,
                self.index,
                seq1.encode(),
                len(seq1),
                seq2.encode(),
                len(seq2),
                avg,
                std,
                imin,
                imax,
            )
        if pe == ffi.NULL:
            return []
        try:
            r1hits = self._conv_hits_raw(pe.regs[0], seq1, min_mapq, min_blen, min_mlen)
            r2hits = self._conv_hits_raw(pe.regs[1], seq2, min_mapq, min_blen, min_mlen)
            results = []
            if r1hits and r2hits:
                c1 = {}
                for h1 in r1hits:
                    if h1[0] not in c1:
                        c1[h1[0]] = []
                    c1[h1[0]].append(h1)
                for h2 in r2hits:
                    if h2[0] in c1:
                        for h1 in c1[h2[0]]:
                            dist = max(h1[2], h2[2]) - min(h1[1], h2[1])
                            if dist < 1000:
                                results.append((h1, h2, True, dist))
            if not results:
                for h1 in r1hits:
                    results.append((h1, None, False, 0))
                for h2 in r2hits:
                    results.append((None, h2, False, 0))
            return results
        finally:
            libbwa.free_pe_regs(pe)

    def _conv_hits_raw(self, regs, seq, min_q, min_blen, min_mlen):
        if regs.n == 0:
            return []
        res_v = libbwa.bwa_mem_reg2aln_all(
            self.opt,
            self.index.bns,
            self.index.pac,
            len(seq),
            seq.encode(),
            ffi.addressof(regs),
            min_q,
            min_blen,
            min_mlen,
        )
        if res_v == ffi.NULL:
            return []
        hits = []
        try:
            for i in range(res_v.n):
                h = res_v.hits[i]
                c_str = ffi.string(h.cigar).decode()
                # (ctg, pos, r_en, strand, q_st, q_en, mapq, cigar, NM, score, rid)
                hits.append(
                    (
                        self._rid_to_name[h.rid],
                        h.pos,
                        h.r_en,
                        h.strand,
                        h.q_st,
                        h.q_en,
                        h.mapq,
                        c_str,
                        h.NM,
                        h.score,
                        h.rid,
                    )
                )
        finally:
            libbwa.free_raw_hit_v(res_v)
        return hits


class BwaIndexer:
    def __init__(
        self, algorithm="auto", block_size=10000000, capture_progress=True, verbose=1
    ):
        self.algorithm, self.block_size, self.capture_progress, self.verbose = (
            algorithm,
            block_size,
            capture_progress,
            verbose,
        )
        try:
            self.algo_type = {"auto": 0, "rb2": 1, "bwtsw": 2, "is": 3}[algorithm]
        except KeyError:
            raise KeyError(f"Unknown algorithm: {algorithm}")
        self.progress = {
            "status": "idle",
            "text_length": 0,
            "iterations": 0,
            "characters_processed": 0,
            "messages": [],
        }

    def build_index(self, fasta_file, prefix=None):
        if not os.path.exists(fasta_file):
            raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
        prefix = prefix or os.path.splitext(fasta_file)[0]
        self.progress.update(
            {
                "status": "building",
                "text_length": 0,
                "iterations": 0,
                "characters_processed": 0,
                "messages": [],
            }
        )
        try:
            libbwa.bwa_verbose = int(self.verbose)
        except:
            pass
        if self.capture_progress:
            rf, wf = os.pipe()
            old_err = os.dup(2)
            os.dup2(wf, 2)

            def reader():
                try:
                    with os.fdopen(rf, "r", errors="replace") as f:
                        for line in f:
                            l = line.rstrip()
                            self.progress["messages"].append(l)
                            m = re.search(r"textLength=(\d+)", l)
                            if m:
                                self.progress["text_length"] = int(m.group(1))
                            m = re.search(
                                r"(\d+) iterations done\. (\d+) characters processed", l
                            )
                            if m:
                                (
                                    self.progress["iterations"],
                                    self.progress["characters_processed"],
                                ) = int(m.group(1)), int(m.group(2))
                except OSError:
                    pass

            t = threading.Thread(target=reader, daemon=True)
            t.start()
            res = libbwa.bwa_idx_build(
                fasta_file.encode(), prefix.encode(), self.algo_type, self.block_size
            )
            # Restore stderr immediately
            os.dup2(old_err, 2)
            os.close(old_err)
            # Close only the write-end to signal EOF to the reader thread
            os.close(wf)
            # Wait for reader to finish and close rf
            t.join(timeout=5)
        else:
            res = libbwa.bwa_idx_build(
                fasta_file.encode(), prefix.encode(), self.algo_type, self.block_size
            )
        if res != 0:
            self.progress["status"] = "failed"
            raise RuntimeError(f"BWA index build failed for {fasta_file}")
        self.progress["status"] = "completed"
        return prefix

    @property
    def progress_percent(self):
        return (
            (self.progress["characters_processed"] / self.progress["text_length"] * 100)
            if self.progress["text_length"] > 0
            else None
        )
