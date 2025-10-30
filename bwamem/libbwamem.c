#include <Python.h>
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "bwamem.h"

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
  int take_all = opt->flag | MEM_F_ALL;
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
  uint8_t* enc = (uint8_t*)malloc(len * sizeof(uint8_t));
  if (enc == NULL) return NULL;
  
  int i;
  for (i = 0; i < len; ++i) {
    // Use BWA's native encoding table
    enc[i] = nst_nt4_table[(int)seq[i]];
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
  
  static const char op_chars[] = "MIDSH";
  char* ptr = result;
  int i;
  
  for (i = 0; i < n_cigar; ++i) {
    uint32_t op_len = cigar[i] >> 4;
    uint32_t op = cigar[i] & 0xF;
    
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
    *ptr++ = op_chars[op];
  }
  
  *ptr = '\0';
  return result;
}
