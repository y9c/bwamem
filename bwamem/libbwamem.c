#include <Python.h>
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "bwamem.h"
#include "../bwa/kstring.h"

// External BWA global variable for verbosity control
extern int bwa_verbose;

// External BWA functions not in bwamem.h
void mem_matesw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], const mem_alnreg_t *a, int l_ms, const uint8_t *ms, mem_alnreg_v *regs);
int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id);
extern unsigned char nst_nt4_table[256];

static PyMethodDef module_functions[] = {{NULL, NULL, 0, NULL}};

#if PY_MAJOR_VERSION >= 3
#define MOD_ERROR_VAL NULL
#define MOD_SUCCESS_VAL(val) val
#define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#define MOD_DEF(ob, name, doc, methods)   \
  static struct PyModuleDef moduledef = { \
      PyModuleDef_HEAD_INIT,              \
      name,                               \
      doc,                                \
      -1,                                 \
      methods,                            \
  };                                      \
  ob = PyModule_Create(&moduledef);
#else
#define MOD_ERROR_VAL
#define MOD_SUCCESS_VAL(val)
#define MOD_INIT(name) void init##name(void)
#define MOD_DEF(ob, name, doc, methods) ob = Py_InitModule3(name, methods, doc);
#endif

MOD_INIT(bwamemy) {
  PyObject* m;
  bwa_verbose = 0; // Silence BWA by default
  MOD_DEF(m, "bwalib", "High-level binding to bwa mem", module_functions)
  if (m == NULL) return MOD_ERROR_VAL;
  return MOD_SUCCESS_VAL(m);
}

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#ifdef MODULE_API_EXPORTS
#define MODULE_API __declspec(dllexport)
#define restrict __restrict
#else
#define MODULE_API __declspec(dllimport)
#endif
#else
#if defined(__GNUC__)
#define MODULE_API __attribute__((visibility("default"))) __attribute__((used))
#else
#define MODULE_API
#endif
#endif

MODULE_API int module_init();

#ifdef __cplusplus
}
#endif

// --- BWA internal CIGAR mapping (from bwamem.h) ---
// 0:M, 1:I, 2:D, 3:S, 4:H
static const char BWA_OP_CHARS[] = "MIDSH";

typedef struct {
  size_t n;
  mem_aln_t* aln;
} mem_aln_v;

mem_aln_v* new_mem_aln_v(size_t n) {
  mem_aln_v* alns = malloc(sizeof(mem_aln_v));
  if (alns == NULL) return NULL;
  alns->aln = malloc(n * sizeof(mem_aln_t));
  if (alns->aln == NULL) { free(alns); return NULL; }
  alns->n = n;
  return alns;
}

void free_mem_aln_v(mem_aln_v* alns) {
  if (alns != NULL) {
    for (size_t i = 0; i < alns->n; ++i) {
        free(alns->aln[i].cigar);
    }
    free(alns->aln);
    free(alns);
  }
}

bwaidx_t* bwa_idx_load_all(const char* hint) {
  return bwa_idx_load(hint, BWA_IDX_ALL);
}

size_t count_primary(mem_alnreg_v* ar) {
  size_t primary = 0;
  for (size_t i = 0; i < ar->n; ++i) { if (ar->a[i].secondary < 0) ++primary; }
  return primary;
}

mem_aln_v* align(mem_opt_t* opt, bwaidx_t* idx, char* seq) {
  size_t seq_len = strlen(seq);
  mem_alnreg_v ar = mem_align1(opt, idx->bwt, idx->bns, idx->pac, seq_len, seq);
  int take_all = opt->flag & MEM_F_ALL;
  size_t n_alns = take_all ? ar.n : count_primary(&ar);
  mem_aln_v* alns = n_alns ? new_mem_aln_v(ar.n) : NULL;
  size_t j = 0;
  for (size_t i = 0; i < ar.n; ++i) {
    if (!take_all && ar.a[i].secondary) continue;
    alns->aln[j++] = mem_reg2aln(opt, idx->bns, idx->pac, seq_len, seq, &ar.a[i]);
  }
  free(ar.a);
  return alns;
}

MODULE_API mem_aln_t* mem_reg2aln_ptr(const mem_opt_t* opt, const bntseq_t* bns, const uint8_t* pac, int l_seq, const char* seq, const mem_alnreg_t* ar) {
  mem_aln_t* out = (mem_aln_t*)malloc(sizeof(mem_aln_t));
  if (out == NULL) return NULL;
  *out = mem_reg2aln(opt, bns, pac, l_seq, seq, ar);
  return out;
}

MODULE_API uint8_t* encode_seq(const char* seq, int len) {
  if (seq == NULL || len < 0) return NULL;
  uint8_t* enc = (uint8_t*)malloc((len > 0 ? (size_t)len : 1) * sizeof(uint8_t));
  if (enc == NULL) return NULL;
  for (int i = 0; i < len; ++i) {
      enc[i] = nst_nt4_table[(unsigned char)seq[i]];
  }
  return enc;
}

typedef struct { uint32_t len; uint32_t op; } cigar_pair_t;

MODULE_API cigar_pair_t* build_cigar_array(const uint32_t* cigar, int n_cigar) {
  if (n_cigar == 0 || cigar == NULL) return NULL;
  cigar_pair_t* result = (cigar_pair_t*)malloc((size_t)n_cigar * sizeof(cigar_pair_t));
  if (result == NULL) return NULL;
  for (int i = 0; i < n_cigar; ++i) { result[i].len = cigar[i] >> 4; result[i].op = cigar[i] & 0xF; }
  return result;
}

MODULE_API int get_cigar_str_len(const uint32_t* cigar, int n_cigar) {
  if (n_cigar == 0 || cigar == NULL) return 0;
  int len = 0;
  for (int i = 0; i < n_cigar; ++i) {
    uint32_t op_len = cigar[i] >> 4;
    if (op_len == 0) len += 2;
    else { uint32_t temp = op_len; while (temp > 0) { len++; temp /= 10; } len++; }
  }
  return len;
}

MODULE_API char* build_cigar_string(const uint32_t* cigar, int n_cigar) {
  if (n_cigar == 0 || cigar == NULL) return NULL;
  kstring_t str = {0, 0, 0};
  for (int i = 0; i < n_cigar; ++i) {
    uint32_t op = cigar[i] & 0xF;
    ksprintf(&str, "%u%c", cigar[i] >> 4, op < 5 ? BWA_OP_CHARS[op] : '?');
  }
  return str.s;
}

// Alias for Python compatibility
MODULE_API char* build_cigar_string_fixed(const uint32_t *cigar, int n_cigar, int l_seq) {
    (void)l_seq; // Silence unused warning
    return build_cigar_string(cigar, n_cigar);
}

MODULE_API int64_t compute_r_en(int64_t r_st, const cigar_pair_t* cigar, int n_cigar) {
  int64_t pos = r_st;
  for (int i = 0; i < n_cigar; ++i) {
      if (cigar[i].op == 0 || cigar[i].op == 2) {
          pos += cigar[i].len;
      }
  }
  return pos;
}

MODULE_API int compute_blen(const cigar_pair_t* cigar, int n_cigar) {
  int length = 0;
  for (int i = 0; i < n_cigar; ++i) {
      if (cigar[i].op <= 2) {
          length += (int)cigar[i].len;
      }
  }
  return length;
}

MODULE_API int compute_mlen(const cigar_pair_t* cigar, int n_cigar) {
  int matches = 0;
  for (int i = 0; i < n_cigar; ++i) {
      if (cigar[i].op == 0) {
          matches += (int)cigar[i].len;
      }
  }
  return matches;
}

// Fast C-based sequence extraction from BWA index
MODULE_API char* bwa_fetch_seq(const bwaidx_t* idx, int rid, int64_t start, int64_t end) {
  if (idx == NULL || rid < 0 || rid >= idx->bns->n_seqs) return NULL;
  int64_t seq_len = idx->bns->anns[rid].len, seq_offset = idx->bns->anns[rid].offset;
  if (start < 0) start = 0;
  if (end > seq_len) end = seq_len;
  if (start >= end) return NULL;
  int64_t len = end - start;
  char* result = (char*)malloc((size_t)len + 1);
  if (result == NULL) return NULL;
  static const char base_chars[] = "ACGT";
  for (int64_t i = 0; i < len; ++i) {
    int64_t gp = seq_offset + start + i;
    result[i] = base_chars[(idx->pac[gp >> 2] >> ((3 - (gp & 3)) << 1)) & 3];
  }
  result[len] = '\0';
  return result;
}

typedef struct { char* cigar; long long pos; long long r_en; int rid; int strand; int q_st; int q_en; int mapq; int NM; int score; } raw_hit_t;
typedef struct { size_t n; raw_hit_t* hits; } raw_hit_v;

MODULE_API void free_raw_hit_v(raw_hit_v* v) {
  if (v) {
      for (size_t i = 0; i < v->n; i++) {
          free(v->hits[i].cigar);
      }
      free(v->hits);
      free(v);
  }
}

MODULE_API raw_hit_v* bwa_mem_reg2aln_all(const mem_opt_t* opt, const bntseq_t* bns, const uint8_t* pac, int l_seq, const char* seq, const mem_alnreg_v* regs, int min_mapq, int min_blen, int min_mlen) {
  if (!regs || regs->n == 0) return NULL;
  raw_hit_v* result = malloc(sizeof(raw_hit_v));
  result->hits = malloc(regs->n * sizeof(raw_hit_t)); result->n = 0;
  for (size_t i = 0; i < regs->n; i++) {
    if (regs->a[i].score < opt->T) continue;
    mem_aln_t aln = mem_reg2aln(opt, bns, pac, l_seq, seq, &regs->a[i]);
    if (aln.rid >= 0 && aln.mapq >= min_mapq) {
      raw_hit_t* h = &result->hits[result->n];
      h->cigar = build_cigar_string(aln.cigar, aln.n_cigar);
      h->pos = aln.pos; h->rid = aln.rid; h->strand = aln.is_rev ? -1 : 1;
      h->q_st = regs->a[i].qb; h->q_en = regs->a[i].qe; h->mapq = aln.mapq; h->NM = aln.NM; h->score = aln.score;
      h->r_en = aln.pos;
      for (int j = 0; j < aln.n_cigar; j++) {
        uint32_t op = aln.cigar[j] & 0xF;
        if (op == 0 || op == 2) h->r_en += (aln.cigar[j] >> 4);
      }
      result->n++;
    }
    free(aln.cigar);
  }
  return result;
}

typedef struct { mem_alnreg_v regs[2]; mem_pestat_t pes[4]; } pe_regs_t;

MODULE_API pe_regs_t* bwa_align_pe(const mem_opt_t* opt, const bwaidx_t* idx, const char* seq1, int l1, const char* seq2, int l2, double avg, double std, int imin, int imax) {
    bwa_verbose = 0;
    pe_regs_t* result = calloc(1, sizeof(pe_regs_t));
    result->regs[0] = mem_align1(opt, idx->bwt, idx->bns, idx->pac, l1, seq1);
    result->regs[1] = mem_align1(opt, idx->bwt, idx->bns, idx->pac, l2, seq2);
    mem_pestat_t* pes = result->pes;
    if (avg > 0) {
        pes[1].failed = 0; pes[1].avg = avg; pes[1].std = std;
        pes[1].low = imin > 0 ? imin : (int)(avg - 4*std);
        pes[1].high = imax > 0 ? imax : (int)(avg + 4*std);
    } else if (result->regs[0].n > 0 && result->regs[1].n > 0) {
        mem_pestat(opt, idx->bns->l_pac, 2, result->regs, pes);
    }
    for (int r = 0; r < 4; r++) {
        if (pes[r].low == 0 && pes[r].high <= 1) {
            pes[r].low = 20; pes[r].high = 600; pes[r].failed = 0; pes[r].avg = 200; pes[r].std = 100;
        }
    }
    if (!(opt->flag & 0x20)) {
        uint8_t *s1e = encode_seq(seq1, l1), *s2e = encode_seq(seq2, l2);
        size_t n1 = (size_t)opt->max_matesw < result->regs[0].n ? (size_t)opt->max_matesw : result->regs[0].n;
        for (size_t j = 0; j < n1; j++) {
            if (result->regs[0].a[j].score >= result->regs[0].a[0].score - opt->pen_unpaired) {
                mem_matesw(opt, idx->bns, idx->pac, pes, &result->regs[0].a[j], l2, s2e, &result->regs[1]);
            }
        }
        size_t n2 = (size_t)opt->max_matesw < result->regs[1].n ? (size_t)opt->max_matesw : result->regs[1].n;
        for (size_t j = 0; j < n2; j++) {
            if (result->regs[1].a[j].score >= result->regs[1].a[0].score - opt->pen_unpaired) {
                mem_matesw(opt, idx->bns, idx->pac, pes, &result->regs[1].a[j], l1, s1e, &result->regs[0]);
            }
        }
        free(s1e); free(s2e);
    }
    mem_mark_primary_se(opt, (int)result->regs[0].n, result->regs[0].a, 0);
    mem_mark_primary_se(opt, (int)result->regs[1].n, result->regs[1].a, 1);
    return result;
}

MODULE_API void free_pe_regs(pe_regs_t* p) {
    if (p) {
        if (p->regs[0].a) free(p->regs[0].a);
        if (p->regs[1].a) free(p->regs[1].a);
        free(p);
    }
}

// --- Original Visualization functions (restored with correct mapping) ---

typedef struct { kstring_t ref_aligned; kstring_t query_aligned; } aligned_seq_pair_t;

MODULE_API void free_aligned_seq_pair(aligned_seq_pair_t* pair) {
    if (pair) {
        free(pair->ref_aligned.s);
        free(pair->query_aligned.s);
        free(pair);
    }
}

MODULE_API aligned_seq_pair_t* apply_cigar_to_sequences(const char* ref_seq, int ref_len, const char* query_seq, int query_len, const cigar_pair_t* cigar, int n_cigar, int q_st, int q_en) {
  if (!ref_seq || !query_seq || !cigar) return NULL;
  aligned_seq_pair_t* res = calloc(1, sizeof(aligned_seq_pair_t));
  int rp = 0, qp = q_st;
  for (int i = 0; i < n_cigar; ++i) {
    uint32_t l = cigar[i].len, op = cigar[i].op;
    if (op == 0) { // M
      int ml = (rp + (int)l > ref_len) ? ref_len - rp : (int)l;
      if (qp + ml > q_en) ml = q_en - qp;
      if (ml > 0) {
          kputsn(ref_seq + rp, ml, &res->ref_aligned);
          kputsn(query_seq + qp, ml, &res->query_aligned);
          rp += (int)l; qp += (int)l;
      }
    } else if (op == 1) { // I
      if (qp < q_en) {
          int il = (qp + (int)l > q_en) ? q_en - qp : (int)l;
          int gl = (int)l;
          while (gl--) kputc('-', &res->ref_aligned);
          kputsn(query_seq + qp, il, &res->query_aligned);
          qp += (int)l;
      }
    } else if (op == 2) { // D
      if (rp < ref_len) {
          int dl = (rp + (int)l > ref_len) ? ref_len - rp : (int)l;
          kputsn(ref_seq + rp, dl, &res->ref_aligned);
          int gl = (int)l;
          while (gl--) kputc('-', &res->query_aligned);
          rp += (int)l;
      }
    } else if (op == 3) {
        qp += (int)l; // S
    }
  }
  return res;
}

MODULE_API char* build_match_indicators(const char* ra, const char* qa, int len) {
  kstring_t s = {0,0,0};
  for (int i = 0; i < len; i++) {
      kputc((ra[i] == '-' || qa[i] == '-') ? ' ' : (ra[i] == qa[i] ? '|' : 'X'), &s);
  }
  return s.s;
}

MODULE_API char* reverse_complement_seq(const char* seq, int len) {
  if (!seq || len <= 0) return NULL;
  char* res = malloc((size_t)len + 1);
  static char ct[256];
  static int init = 0;
  if (!init) {
      for (int i = 0; i < 256; i++) ct[i] = (char)i;
      ct['A'] = 'T'; ct['C'] = 'G'; ct['G'] = 'C'; ct['T'] = 'A';
      ct['a'] = 't'; ct['c'] = 'g'; ct['g'] = 'c'; ct['t'] = 'a';
      init = 1;
  }
  for (int i = 0; i < len; i++) {
      res[i] = ct[(unsigned char)seq[len - 1 - i]];
  }
  res[len] = '\0';
  return res;
}

MODULE_API char* visualize_alignment_c(const char* ctg, int64_t r_st, int64_t r_en, int strand, int score, int mapq, const char* ref_seq, int ref_len, const char* query_seq, int query_len, const cigar_pair_t* cigar, int n_cigar, int q_st, int q_en, int lw) {
  char* qsw = (strand < 0) ? reverse_complement_seq(query_seq, query_len) : NULL;
  if (qsw) query_seq = qsw;
  aligned_seq_pair_t* ap = apply_cigar_to_sequences(ref_seq, ref_len, query_seq, query_len, cigar, n_cigar, q_st, q_en);
  if (!ap) { if(qsw) free(qsw); return NULL; }
  char* ml = build_match_indicators(ap->ref_aligned.s, ap->query_aligned.s, (int)ap->ref_aligned.l);
  kstring_t out = {0,0,0};
  ksprintf(&out, "%s:%ld-%ld (strand: %s, score: %d, mapq: %d)\n", ctg, (long)r_st, (long)r_en, strand>0?"+":"-", score, mapq);
  int al = (int)ap->ref_aligned.l, nc = (al+lw-1)/lw;
  for (int c=0; c<nc; c++) {
    int cs = c*lw, cl = (cs+lw < al) ? lw : al-cs;
    int64_t rs = r_st, qs = q_st;
    for (int j=0; j<cs; j++) { if(ap->ref_aligned.s[j]!='-') rs++; if(ap->query_aligned.s[j]!='-') qs++; }
    int64_t re = rs, qe = qs;
    for (int j=cs; j<cs+cl; j++) { if(ap->ref_aligned.s[j]!='-') re++; if(ap->query_aligned.s[j]!='-') qe++; }
    ksprintf(&out, "Ref %6ld ", (long)rs); kputsn(ap->ref_aligned.s+cs, cl, &out); ksprintf(&out, " %ld\n", (long)re-1);
    kputs("            ", &out); kputsn(ml+cs, cl, &out); kputc('\n', &out);
    ksprintf(&out, "Qry %6d ", (int)qs); kputsn(ap->query_aligned.s+cs, cl, &out); ksprintf(&out, " %d\n", (int)qe-1);
    if (c<nc-1) kputc('\n', &out);
  }
  free(ml); free_aligned_seq_pair(ap); if(qsw) free(qsw); return out.s;
}
