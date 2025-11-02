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

  // Set BWA verbosity to 1 (quiet mode - only warnings/errors)
  bwa_verbose = 1;

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
  
  aligned_seq_pair_t* result = (aligned_seq_pair_t*)calloc(1, sizeof(aligned_seq_pair_t));
  if (result == NULL) return NULL;
  
  int ref_pos = 0;
  int query_pos = q_st;
  
  int i;
  for (i = 0; i < n_cigar; ++i) {
    uint32_t length = cigar[i].len;
    uint32_t op = cigar[i].op;
    
    switch (op) {
      case 0:  // M (match/mismatch)
        if (ref_pos + length <= ref_len && query_pos + length <= query_len) {
          kputsn(ref_seq + ref_pos, length, &result->ref_aligned);
          kputsn(query_seq + query_pos, length, &result->query_aligned);
          ref_pos += length;
          query_pos += length;
        }
        break;
        
      case 1:  // I (insertion in query)
        if (query_pos + length <= query_len) {
          int gap_len = length;
          while (gap_len-- > 0) kputc('-', &result->ref_aligned);
          kputsn(query_seq + query_pos, length, &result->query_aligned);
          query_pos += length;
        }
        break;
        
      case 2:  // D (deletion in query)
        if (ref_pos + length <= ref_len) {
          kputsn(ref_seq + ref_pos, length, &result->ref_aligned);
          int gap_len = length;
          while (gap_len-- > 0) kputc('-', &result->query_aligned);
          ref_pos += length;
        }
        break;
        
      case 4:  // S (soft clip) - skip in query
        query_pos += length;
        break;
        
      case 5:  // H (hard clip) - skip
        break;
        
      default:
        break;
    }
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
  
  // Step 3: Build match indicators
  char* match_line = build_match_indicators(
      aligned_pair->ref_aligned.s,
      aligned_pair->query_aligned.s,
      aligned_pair->ref_aligned.l);
  
  if (match_line == NULL) {
    free_aligned_seq_pair(aligned_pair);
    if (query_seq_work) free(query_seq_work);
    return NULL;
  }
  
  // Step 4: Build output using kstring_t
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
