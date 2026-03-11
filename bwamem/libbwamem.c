#include <Python.h>
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "bwamem.h"
#include "../bwa/kstring.h"

// External BWA functions not in bwamem.h
void mem_matesw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], const mem_alnreg_t *a, int l_ms, const uint8_t *ms, mem_alnreg_v *regs);
int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id);

// External BWA global variable for verbosity control
extern int bwa_verbose;

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

  // Set BWA verbosity to 0 (silent mode) by default
  // This can be adjusted by users if needed, but silent is better for library usage
  bwa_verbose = 0;

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

typedef struct {
  size_t n;
  mem_aln_t* aln;
} mem_aln_v;

mem_aln_v* new_mem_aln_v(size_t n) {
  // Allocate a mem_aln_t vector
  mem_aln_v* alns = malloc(sizeof(mem_aln_v));
  if (alns == NULL) return NULL;
  alns->aln = malloc(n * sizeof(mem_aln_t));
  if (alns->aln == NULL) {
    free(alns);
    return NULL;
  }
  alns->n = n;
  return alns;
}

void free_mem_aln_v(mem_aln_v* alns) {
  // free mem_aln_v and all its submembers
  if (alns != NULL) {
    for (size_t i = 0; i < alns->n; ++i) {
      free(alns->aln[i].cigar);
    }
    free(alns->aln);
    free(alns);
  }
}

// (no wrapper needed; link whole lib to export bwa_idx_build)

bwaidx_t* bwa_idx_load_all(const char* hint) {
  return bwa_idx_load(hint, BWA_IDX_ALL);
}

size_t count_primary(mem_alnreg_v* ar) {
  size_t primary = 0;
  for (size_t i = 0; i < ar->n; ++i) {
    if (ar->a[i].secondary >= 0) continue;
    ++primary;
  }
  return primary;
}

mem_aln_v* align(mem_opt_t* opt, bwaidx_t* idx, char* seq) {
  // call the aligner
  size_t seq_len = strlen(seq);
  mem_alnreg_v ar = mem_align1(opt, idx->bwt, idx->bns, idx->pac, seq_len, seq);

  // check if we take all or only primary alignments
  int take_all = opt->flag & MEM_F_ALL;
  size_t n_alns = take_all ? ar.n : count_primary(&ar);

  // allocate memory for the result if there any
  mem_aln_v* alns = n_alns ? new_mem_aln_v(ar.n) : NULL;

  // copy results (if there are any)
  size_t j = 0;
  for (size_t i = 0; i < ar.n; ++i) {
    if (!take_all && ar.a[i].secondary) continue;
    alns->aln[j++] =
        mem_reg2aln(opt, idx->bns, idx->pac, seq_len, seq, &ar.a[i]);
  }

  // free the intermidiate results and return
  free(ar.a);
  return alns;
}

// Wrapper to return a pointer to mem_aln_t to avoid returning
// a struct-with-bitfields directly across the FFI boundary
MODULE_API mem_aln_t* mem_reg2aln_ptr(
    const mem_opt_t* opt,
    const bntseq_t* bns,
    const uint8_t* pac,
    int l_seq,
    const char* seq,
    const mem_alnreg_t* ar) {
  mem_aln_t* out = (mem_aln_t*)malloc(sizeof(mem_aln_t));
  if (out == NULL) return NULL;
  *out = mem_reg2aln(opt, bns, pac, l_seq, seq, ar);
  return out;
}

// BWA's nucleotide encoding table
extern unsigned char nst_nt4_table[256];

// Fast C-based sequence encoding using BWA's native table
// Converts ASCII sequence to 0-3 encoding (A=0, C=1, G=2, T=3, N=4)
// Returns allocated uint8_t array that caller must free()
MODULE_API uint8_t* encode_seq(const char* seq, int len) {
  // Input validation
  if (seq == NULL || len < 0) return NULL;
  if (len == 0) {
    // Return empty allocation for zero-length sequences
    uint8_t* enc = (uint8_t*)malloc(1);
    return enc;
  }
  
  uint8_t* enc = (uint8_t*)malloc(len * sizeof(uint8_t));
  if (enc == NULL) return NULL;
  
  int i;
  for (i = 0; i < len; ++i) {
    // Use BWA's native encoding table
    // The table safely handles all 256 possible byte values
    enc[i] = nst_nt4_table[(unsigned char)seq[i]];
  }
  
  return enc;
}

// Structure to hold CIGAR operation pairs [length, op]
typedef struct {
  uint32_t len;
  uint32_t op;
} cigar_pair_t;

// Fast C-based CIGAR array builder
// Converts BAM-encoded CIGAR (opLen<<4|op) to array of [length, op] pairs
// Returns allocated array that caller must free()
MODULE_API cigar_pair_t* build_cigar_array(const uint32_t* cigar, int n_cigar) {
  if (n_cigar == 0 || cigar == NULL) return NULL;
  
  cigar_pair_t* result = (cigar_pair_t*)malloc(n_cigar * sizeof(cigar_pair_t));
  if (result == NULL) return NULL;
  
  int i;
  for (i = 0; i < n_cigar; ++i) {
    result[i].len = cigar[i] >> 4;
    result[i].op = cigar[i] & 0xF;
  }
  
  return result;
}

// Get CIGAR string length (for pre-allocation)
MODULE_API int get_cigar_str_len(const uint32_t* cigar, int n_cigar) {
  if (n_cigar == 0 || cigar == NULL) return 0;
  
  int len = 0;
  int i;
  for (i = 0; i < n_cigar; ++i) {
    uint32_t op_len = cigar[i] >> 4;
    // Count digits in op_len + 1 for operation character
    if (op_len == 0) {
      len += 2; // "0" + op
    } else {
      uint32_t temp = op_len;
      while (temp > 0) {
        len++;
        temp /= 10;
      }
      len++; // for operation character
    }
  }
  return len;
}

// Build CIGAR string directly in C (faster than Python string building)
MODULE_API char* build_cigar_string(const uint32_t* cigar, int n_cigar) {
  if (n_cigar == 0 || cigar == NULL) return NULL;
  
  // Pre-calculate length needed
  int len = get_cigar_str_len(cigar, n_cigar);
  char* result = (char*)malloc((len + 1) * sizeof(char)); // +1 for null terminator
  if (result == NULL) return NULL;
  
  // Extended op_chars to handle all possible CIGAR operations
  // BAM spec: M=0, I=1, D=2, N=3, S=4, H=5, P=6, =7, X=8, B=9
  static const char op_chars[] = "MIDNSHP=XB";
  char* ptr = result;
  int i;
  
  for (i = 0; i < n_cigar; ++i) {
    uint32_t op_len = cigar[i] >> 4;
    uint32_t op = cigar[i] & 0xF;
    
    // Bounds check for operation code
    if (op >= sizeof(op_chars) - 1) {
      // Invalid operation code - use '?' as placeholder
      op = '?';
    }
    
    // Convert length to string
    if (op_len == 0) {
      *ptr++ = '0';
    } else {
      // Convert number to string in reverse
      char temp[20];
      int pos = 0;
      uint32_t temp_len = op_len;
      while (temp_len > 0) {
        temp[pos++] = '0' + (temp_len % 10);
        temp_len /= 10;
      }
      // Copy reversed digits
      while (pos > 0) {
        *ptr++ = temp[--pos];
      }
    }
    
    // Add operation character
    if (op == '?') {
      *ptr++ = '?';
    } else {
      *ptr++ = op_chars[op];
    }
  }
  
  *ptr = '\0';
  return result;
}

// Fast C-based calculation of reference end position from CIGAR
// r_en = r_st + sum of lengths of operations that consume reference (M, D, N)
MODULE_API int64_t compute_r_en(int64_t r_st, const cigar_pair_t* cigar, int n_cigar) {
  if (cigar == NULL || n_cigar == 0) return r_st;
  
  int64_t pos = r_st;
  int i;
  for (i = 0; i < n_cigar; ++i) {
    uint32_t op = cigar[i].op;
    // M=0, D=2, N=3 consume reference
    if (op == 0 || op == 2 || op == 3) {
      pos += cigar[i].len;
    }
  }
  return pos;
}

// Fast C-based calculation of alignment block length (including gaps)
// blen = sum of lengths of M, I, D, N operations
MODULE_API int compute_blen(const cigar_pair_t* cigar, int n_cigar) {
  if (cigar == NULL || n_cigar == 0) return 0;
  
  int length = 0;
  int i;
  for (i = 0; i < n_cigar; ++i) {
    uint32_t op = cigar[i].op;
    // M=0, I=1, D=2, N=3
    if (op <= 3) {
      length += cigar[i].len;
    }
  }
  return length;
}

// Fast C-based calculation of number of matching bases
// mlen = sum of lengths of M operations
MODULE_API int compute_mlen(const cigar_pair_t* cigar, int n_cigar) {
  if (cigar == NULL || n_cigar == 0) return 0;
  
  int matches = 0;
  int i;
  for (i = 0; i < n_cigar; ++i) {
    // M=0 (match/mismatch)
    if (cigar[i].op == 0) {
      matches += cigar[i].len;
    }
  }
  return matches;
}

// Fast C-based sequence extraction from BWA index
// Returns allocated string that caller must free()
MODULE_API char* bwa_fetch_seq(const bwaidx_t* idx, int rid, int64_t start, int64_t end) {
  if (idx == NULL || rid < 0 || rid >= idx->bns->n_seqs) return NULL;

  int64_t seq_len = idx->bns->anns[rid].len;
  int64_t seq_offset = idx->bns->anns[rid].offset;

  if (start < 0) start = 0;
  if (end > seq_len) end = seq_len;
  if (start >= end) return NULL;

  
  int64_t len = end - start;
  char* result = (char*)malloc(len + 1);
  if (result == NULL) return NULL;
  
  static const char base_chars[] = "ACGT";
  int64_t i;
  for (i = 0; i < len; ++i) {
    int64_t global_pos = seq_offset + start + i;
    int64_t byte_pos = global_pos >> 2;
    int bit_offset = (3 - (global_pos & 3)) << 1;
    result[i] = base_chars[(idx->pac[byte_pos] >> bit_offset) & 3];
  }
  result[len] = '\0';
  return result;
}

// Structure to hold aligned sequences for visualization
// Using kstring_t for automatic memory management
typedef struct {
  kstring_t ref_aligned;    // Reference sequence with gaps
  kstring_t query_aligned;  // Query sequence with gaps
} aligned_seq_pair_t;

// Free aligned sequence pair
MODULE_API void free_aligned_seq_pair(aligned_seq_pair_t* pair) {
  if (pair != NULL) {
    free(pair->ref_aligned.s);
    free(pair->query_aligned.s);
    free(pair);
  }
}

// Reverse complement a DNA sequence in C
MODULE_API char* reverse_complement_seq(const char* seq, int len) {
  if (seq == NULL || len <= 0) return NULL;
  
  char* result = (char*)malloc((len + 1) * sizeof(char));
  if (result == NULL) return NULL;
  
  // Complement table - initialized to 0, then set specific values
  static char comp_table[256] = {0};
  static int table_initialized = 0;
  
  if (!table_initialized) {
    comp_table['A'] = 'T'; comp_table['a'] = 't';
    comp_table['T'] = 'A'; comp_table['t'] = 'a';
    comp_table['G'] = 'C'; comp_table['g'] = 'c';
    comp_table['C'] = 'G'; comp_table['c'] = 'g';
    comp_table['N'] = 'N'; comp_table['n'] = 'n';
    table_initialized = 1;
  }
  
  int i;
  for (i = 0; i < len; ++i) {
    char base = seq[len - 1 - i];  // Reverse
    unsigned char idx = (unsigned char)base;
    result[i] = comp_table[idx] ? comp_table[idx] : base;
  }
  result[len] = '\0';
  return result;
}

// Apply CIGAR operations to create aligned sequences (in C for performance)
// Returns aligned sequence pair or NULL on error
MODULE_API aligned_seq_pair_t* apply_cigar_to_sequences(
    const char* ref_seq, int ref_len,
    const char* query_seq, int query_len,
    const cigar_pair_t* cigar, int n_cigar,
    int q_st, int q_en) {
  
  if (ref_seq == NULL || query_seq == NULL || cigar == NULL || n_cigar == 0) {
    return NULL;
  }
  
  if (q_st < 0 || q_en < 0 || q_st >= query_len || q_en > query_len || q_st > q_en) {
    return NULL;
  }
  
  aligned_seq_pair_t* result = (aligned_seq_pair_t*)calloc(1, sizeof(aligned_seq_pair_t));
  if (result == NULL) return NULL;
  
  // Initialize kstring_t structures
  result->ref_aligned = (kstring_t){0, 0, NULL};
  result->query_aligned = (kstring_t){0, 0, NULL};
  
  int ref_pos = 0;
  int query_pos = q_st;
  
  int i;
  for (i = 0; i < n_cigar; ++i) {
    uint32_t length = cigar[i].len;
    uint32_t op = cigar[i].op;
    
    switch (op) {
      case 0:  // M (match/mismatch)
        // Check bounds - use q_en as the limit for query sequence
        if (ref_pos < ref_len && query_pos < q_en) {
          int match_len = length;
          if (ref_pos + match_len > ref_len) match_len = ref_len - ref_pos;
          if (query_pos + match_len > q_en) match_len = q_en - query_pos;
          if (match_len > 0) {
            kputsn(ref_seq + ref_pos, match_len, &result->ref_aligned);
            kputsn(query_seq + query_pos, match_len, &result->query_aligned);
            ref_pos += length;  // Use original length for ref_pos tracking
            query_pos += length;  // Use original length for query_pos tracking
          }
        }
        break;
        
      case 1:  // I (insertion in query)
        if (query_pos < q_en) {
          int ins_len = length;
          if (query_pos + ins_len > q_en) ins_len = q_en - query_pos;
          if (ins_len > 0) {
            int gap_len = length;  // Always use full length for gaps
            while (gap_len-- > 0) kputc('-', &result->ref_aligned);
            kputsn(query_seq + query_pos, ins_len, &result->query_aligned);
            query_pos += length;  // Use original length
          }
        }
        break;
        
      case 2:  // D (deletion in query)
        if (ref_pos < ref_len) {
          int del_len = length;
          if (ref_pos + del_len > ref_len) del_len = ref_len - ref_pos;
          if (del_len > 0) {
            kputsn(ref_seq + ref_pos, del_len, &result->ref_aligned);
            int gap_len = length;  // Always use full length for gaps
            while (gap_len-- > 0) kputc('-', &result->query_aligned);
            ref_pos += length;  // Use original length
          }
        }
        break;
        
      case 4:  // S (soft clip) - skip in query (don't add to alignment)
        // Soft clip doesn't consume reference, just advances query position
        query_pos += length;
        break;
        
      case 5:  // H (hard clip) - skip
        break;
        
      default:
        break;
    }
  }
  
  // Validate that we actually produced some aligned sequences
  if (result->ref_aligned.l == 0 || result->query_aligned.l == 0) {
    free_aligned_seq_pair(result);
    return NULL;
  }
  
  return result;
}

// Build match indicator line from aligned sequences (| for match, X for mismatch, space for gap)
MODULE_API char* build_match_indicators(const char* ref_aligned, const char* query_aligned, int len) {
  if (ref_aligned == NULL || query_aligned == NULL || len <= 0) return NULL;
  
  kstring_t str = {0, 0, NULL};
  int i;
  for (i = 0; i < len; ++i) {
    char r = ref_aligned[i];
    char q = query_aligned[i];
    kputc((r == '-' || q == '-') ? ' ' : (r == q ? '|' : 'X'), &str);
  }
  return str.s;  // Caller must free
}


// Structure for directional mapping results
typedef struct {
  char* cigar;
  char* md;
  long long pos;
  int rid;
  int strand;
  int score;
  int bad_mm;
  int yf, zf, yc, zc, ns, nc;
} dir_hit_t;

typedef struct {
  size_t n;
  dir_hit_t* hits;
} dir_hit_v;

MODULE_API void free_dir_hit_v(dir_hit_v* v) {
  if (v) {
    for (size_t i = 0; i < v->n; i++) {
      free(v->hits[i].cigar);
      free(v->hits[i].md);
    }
    free(v->hits);
    free(v);
  }
}

// Internal helper for directional scoring inside the mapping loop
static void score_hit_directional(const bwaidx_t* idx, const char* q_seq, int l_seq, 
                                 const uint32_t* cigar, int n_cigar, int rid, int64_t pos, int is_rev,
                                 int is_orientation1, dir_hit_t* out) {
  // Build CIGAR string
  out->cigar = build_cigar_string(cigar, n_cigar);
  
  // Get reference end pos to know fetch length
  int64_t r_en = pos;
  for (int i = 0; i < n_cigar; i++) {
    uint32_t op = cigar[i] & 0xF;
    if (op == 0 || op == 2 || op == 3) r_en += (cigar[i] >> 4);
  }
  
  // Fetch reference
  char* ref = bwa_fetch_seq(idx, rid, pos, r_en);
  if (!ref) {
    out->score = -1000;
    out->md = strdup("0");
    return;
  }
  
  // Directional scoring logic
  size_t buffer_size = l_seq * 16 + 128;
  char* md_buf = malloc(buffer_size);
  
  int yf = 0, zf = 0, yc = 0, zc = 0, ns = 0, nc = 0;
  int matches = 0, expected_conversions = 0, wrong_conversions = 0, other_mismatches = 0, indels = 0;
  int r_idx = 0, q_idx = 0, match_count = 0;
  size_t md_pos = 0;
  
  char b1, b3;
  if (is_orientation1) { b1 = 'A'; b3 = 'C'; } 
  else { b1 = 'T'; b3 = 'G'; }
  
  for (int i = 0; i < n_cigar; i++) {
    uint32_t len = cigar[i] >> 4;
    uint32_t op = cigar[i] & 0xF;
    
    if (op == 0) { // M
      for (uint32_t j = 0; j < len; j++) {
        char rb = ref[r_idx], qb = q_seq[q_idx];
        if (rb == qb) {
          matches++; match_count++;
          if (qb == b1) zf++; else if (qb == b3) zc++;
        } else {
          int n = snprintf(md_buf + md_pos, buffer_size - md_pos, "%d%c", match_count, rb);
          if (n > 0) md_pos += n;
          match_count = 0;
          if (is_orientation1) {
            if ((rb == 'C' && qb == 'T') || (rb == 'A' && qb == 'G')) {
              expected_conversions++; if (qb == 'G') yf++; else yc++;
            } else { wrong_conversions++; ns++; }
          } else {
            if ((rb == 'G' && qb == 'A') || (rb == 'T' && qb == 'C')) {
              expected_conversions++; if (qb == 'A') yf++; else yc++;
            } else { wrong_conversions++; ns++; }
          }
        }
        r_idx++; q_idx++;
      }
    } else if (op == 1 || op == 4) { q_idx += len; nc += len; if (op == 1) indels += len; }
    else if (op == 2) {
      int n = snprintf(md_buf + md_pos, buffer_size - md_pos, "%d^", match_count);
      if (n > 0) md_pos += n;
      for (uint32_t j = 0; j < len; j++) {
        if (md_pos < buffer_size - 1) md_buf[md_pos++] = ref[r_idx++];
        else r_idx++;
      }
      match_count = 0; nc += len; indels += len;
    } else if (op == 3) { r_idx += len; }
  }
  int final_n = snprintf(md_buf + md_pos, buffer_size - md_pos, "%d", match_count);
  if (final_n > 0) md_pos += final_n;
  md_buf[md_pos] = '\0';
  
  out->md = strdup(md_buf);
  free(md_buf);
  
  out->score = matches + expected_conversions - wrong_conversions - other_mismatches - indels;
  out->bad_mm = wrong_conversions + other_mismatches;
  out->yf = yf; out->zf = zf; out->yc = yc; out->zc = zc; out->ns = ns; out->nc = nc;
  out->rid = rid; out->pos = pos; out->strand = is_rev ? -1 : 1;
  
  free(ref);
}

MODULE_API dir_hit_v* bwa_map_directional(const mem_opt_t* opt, const bwaidx_t* idx, 
                                         const char* q_conv, const char* q_orig, 
                                         int is_orientation1) {
  int l_seq = strlen(q_conv);
  mem_alnreg_v ar = mem_align1(opt, idx->bwt, idx->bns, idx->pac, l_seq, q_conv);
  
  if (ar.n == 0) { free(ar.a); return NULL; }
  
  dir_hit_v* result = malloc(sizeof(dir_hit_v));
  result->hits = malloc(ar.n * sizeof(dir_hit_t));
  result->n = 0;
  
  for (size_t i = 0; i < ar.n; i++) {
    mem_aln_t aln = mem_reg2aln(opt, idx->bns, idx->pac, l_seq, q_conv, &ar.a[i]);
    if (aln.rid >= 0) {
      score_hit_directional(idx, q_orig, l_seq, aln.cigar, aln.n_cigar, aln.rid, aln.pos, aln.is_rev, 
                            is_orientation1, &result->hits[result->n]);
      result->n++;
    }
    free(aln.cigar);
  }
  
  free(ar.a);
  return result;
}

// Structure for raw alignment hits
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

MODULE_API void free_raw_hit_v(raw_hit_v* v) {
  if (v) {
    for (size_t i = 0; i < v->n; i++) {
      free(v->hits[i].cigar);
    }
    free(v->hits);
    free(v);
  }
}

typedef struct {
  mem_alnreg_v regs[2];
  mem_pestat_t pes[4];
} pe_regs_t;

MODULE_API pe_regs_t* bwa_align_pe(const mem_opt_t* opt, const bwaidx_t* idx, const char* seq1, int l1, const char* seq2, int l2, double avg, double std, int imin, int imax) {
    pe_regs_t* result = calloc(1, sizeof(pe_regs_t));
    result->regs[0] = mem_align1(opt, idx->bwt, idx->bns, idx->pac, l1, seq1);
    result->regs[1] = mem_align1(opt, idx->bwt, idx->bns, idx->pac, l2, seq2);
    
    mem_pestat_t* pes = result->pes;
    if (avg > 0) {
        pes[1].failed = 0; pes[1].avg = avg; pes[1].std = std;
        pes[1].low = imin > 0 ? imin : (int)(avg - 4*std);
        pes[1].high = imax > 0 ? imax : (int)(avg + 4*std);
    } else if (result->regs[0].n > 0 && result->regs[1].n > 0) {
        // We assume stderr is suppressed externally if needed
        mem_pestat(opt, idx->bns->l_pac, 2, result->regs, pes);
    }
    
    for (int r = 0; r < 4; r++) {
        if (pes[r].low == 0 && pes[r].high <= 1) {
            pes[r].low = 20; pes[r].high = 600; pes[r].failed = 0; pes[r].avg = 200; pes[r].std = 100;
        }
    }
    
    if (!(opt->flag & 0x20)) {
        uint8_t* s1e = encode_seq(seq1, l1);
        uint8_t* s2e = encode_seq(seq2, l2);
        
        int n1 = opt->max_matesw < result->regs[0].n ? opt->max_matesw : result->regs[0].n;
        for (int j = 0; j < n1; j++) {
            if (result->regs[0].a[j].score >= result->regs[0].a[0].score - opt->pen_unpaired) {
                mem_matesw(opt, idx->bns, idx->pac, pes, &result->regs[0].a[j], l2, s2e, &result->regs[1]);
            }
        }
        
        int n2 = opt->max_matesw < result->regs[1].n ? opt->max_matesw : result->regs[1].n;
        for (int j = 0; j < n2; j++) {
            if (result->regs[1].a[j].score >= result->regs[1].a[0].score - opt->pen_unpaired) {
                mem_matesw(opt, idx->bns, idx->pac, pes, &result->regs[1].a[j], l1, s1e, &result->regs[0]);
            }
        }
        free(s1e); free(s2e);
    }
    
    mem_mark_primary_se(opt, result->regs[0].n, result->regs[0].a, 0);
    mem_mark_primary_se(opt, result->regs[1].n, result->regs[1].a, 1);
    
    return result;
}

MODULE_API void free_pe_regs(pe_regs_t* p) {
    if (p) {
        if (p->regs[0].a) free(p->regs[0].a);
        if (p->regs[1].a) free(p->regs[1].a);
        free(p);
    }
}

MODULE_API raw_hit_v* bwa_mem_reg2aln_all(const mem_opt_t* opt, const bntseq_t* bns, const uint8_t* pac, int l_seq, const char* seq, const mem_alnreg_v* regs, int min_mapq, int min_blen, int min_mlen) {
  if (!regs || regs->n == 0) return NULL;
  
  raw_hit_v* result = malloc(sizeof(raw_hit_v));
  if (!result) return NULL;
  result->hits = malloc(regs->n * sizeof(raw_hit_t));
  if (!result->hits) { free(result); return NULL; }
  result->n = 0;
  
  for (size_t i = 0; i < regs->n; i++) {
    if (regs->a[i].score < opt->T) continue;
    
    mem_aln_t aln = mem_reg2aln(opt, bns, pac, l_seq, seq, &regs->a[i]);
    if (aln.rid >= 0 && aln.mapq >= min_mapq) {
      cigar_pair_t* cp = (cigar_pair_t*)aln.cigar;
      
      int skip = 0;
      if (min_blen > 0 && compute_blen(cp, aln.n_cigar) < min_blen) skip = 1;
      if (!skip && min_mlen > 0 && compute_mlen(cp, aln.n_cigar) < min_mlen) skip = 1;
      
      if (!skip) {
        raw_hit_t* h = &result->hits[result->n];
        h->cigar = build_cigar_string(aln.cigar, aln.n_cigar);
        h->pos = aln.pos;
        h->rid = aln.rid;
        h->strand = aln.is_rev ? -1 : 1;
        h->q_st = regs->a[i].qb;
        h->q_en = regs->a[i].qe;
        h->mapq = aln.mapq;
        h->NM = aln.NM;
        h->score = aln.score;
        
        h->r_en = aln.pos;
        for (int j = 0; j < aln.n_cigar; j++) {
          uint32_t op = aln.cigar[j] & 0xF;
          if (op == 0 || op == 2 || op == 3 || op == 7 || op == 8) {
            h->r_en += (aln.cigar[j] >> 4);
          }
        }
        result->n++;
      }
    }
    free(aln.cigar);
  }
  return result;
}

// Complete visualization in C - returns formatted string
MODULE_API char* visualize_alignment_c(
    const char* ctg, int64_t r_st, int64_t r_en, int strand, int score, int mapq,
    const char* ref_seq, int ref_seq_len,
    const char* query_seq, int query_seq_len,
    const cigar_pair_t* cigar, int n_cigar,
    int q_st, int q_en,
    int line_width) {
  
  if (ctg == NULL || ref_seq == NULL || query_seq == NULL || cigar == NULL || n_cigar == 0) {
    return NULL;
  }
  
  if (line_width <= 0) line_width = 80;
  
  // Step 1: Reverse complement query if needed
  char* query_seq_work = NULL;
  if (strand < 0) {
    query_seq_work = reverse_complement_seq(query_seq, query_seq_len);
    if (query_seq_work == NULL) return NULL;
    query_seq = query_seq_work;
  }
  
  // Step 2: Apply CIGAR
  aligned_seq_pair_t* aligned_pair = apply_cigar_to_sequences(
      ref_seq, ref_seq_len,
      query_seq, query_seq_len,
      cigar, n_cigar,
      q_st, q_en);
  
  if (aligned_pair == NULL) {
    if (query_seq_work) free(query_seq_work);
    return NULL;
  }
  
  // Step 3: Validate aligned sequences
  if (aligned_pair->ref_aligned.s == NULL || aligned_pair->query_aligned.s == NULL ||
      aligned_pair->ref_aligned.l == 0 || aligned_pair->query_aligned.l == 0 ||
      aligned_pair->ref_aligned.l != aligned_pair->query_aligned.l) {
    free_aligned_seq_pair(aligned_pair);
    if (query_seq_work) free(query_seq_work);
    return NULL;
  }
  
  // Step 4: Build match indicators
  char* match_line = build_match_indicators(
      aligned_pair->ref_aligned.s,
      aligned_pair->query_aligned.s,
      aligned_pair->ref_aligned.l);
  
  if (match_line == NULL) {
    free_aligned_seq_pair(aligned_pair);
    if (query_seq_work) free(query_seq_work);
    return NULL;
  }
  
  // Step 5: Build output using kstring_t
  kstring_t output = {0, 0, NULL};
  int aln_len = aligned_pair->ref_aligned.l;
  int num_chunks = (aln_len + line_width - 1) / line_width;
  
  // Build header
  const char* strand_str = (strand > 0) ? "+" : "-";
  ksprintf(&output, "%s:%ld-%ld (strand: %s, score: %d, mapq: %d)\n",
           ctg, (long)r_st, (long)r_en, strand_str, score, mapq);
  
  // Build separator line
  int sep_len = output.l - 1;  // Exclude newline
  if (sep_len > line_width) sep_len = line_width;
  int i;
  for (i = 0; i < sep_len; ++i) kputc('=', &output);
  kputc('\n', &output);
  
  // Build alignment lines in chunks
  int chunk;
  for (chunk = 0; chunk < num_chunks; ++chunk) {
    int chunk_start = chunk * line_width;
    int chunk_len = (chunk_start + line_width < aln_len) 
                     ? line_width 
                     : aln_len - chunk_start;
    
    // Calculate positions for this chunk
    int64_t r_pos = r_st;
    int q_pos = q_st;
    
    // Count non-gap bases up to chunk_start
    int j;
    for (j = 0; j < chunk_start && j < aln_len; ++j) {
      if (aligned_pair->ref_aligned.s[j] != '-') ++r_pos;
      if (aligned_pair->query_aligned.s[j] != '-') ++q_pos;
    }
    
    int64_t r_end_pos = r_pos;
    int q_end_pos = q_pos;
    
    // Count non-gap bases in chunk
    for (j = chunk_start; j < chunk_start + chunk_len && j < aln_len; ++j) {
      if (aligned_pair->ref_aligned.s[j] != '-') ++r_end_pos;
      if (aligned_pair->query_aligned.s[j] != '-') ++q_end_pos;
    }
    
    // Reference line: "Ref  r_pos seq r_end_pos\n"
    kputs("Ref  ", &output);
    ksprintf(&output, "%6ld ", (long)r_pos);
    kputsn(aligned_pair->ref_aligned.s + chunk_start, chunk_len, &output);
    ksprintf(&output, " %ld\n", (long)(r_end_pos - 1));
    
    // Match line: pad to align with sequence start (12 chars: "Ref  " + "r_pos ")
    kputs("            ", &output);  // 12 spaces
    kputsn(match_line + chunk_start, chunk_len, &output);
    kputc('\n', &output);
    
    // Query line: "Qry  q_pos seq q_end_pos\n"
    kputs("Qry  ", &output);
    ksprintf(&output, "%6d ", q_pos);
    kputsn(aligned_pair->query_aligned.s + chunk_start, chunk_len, &output);
    ksprintf(&output, " %d\n", q_end_pos - 1);
    
    // Blank line between chunks (except last)
    if (chunk < num_chunks - 1) kputc('\n', &output);
  }
  
  // Cleanup
  free(match_line);
  free_aligned_seq_pair(aligned_pair);
  if (query_seq_work) free(query_seq_work);
  
  return output.s;  // Caller must free
}
