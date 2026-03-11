#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "bwa.h"
#include "bwamem.h"

static void update_a(mem_opt_t* opt, const mem_opt_t* opt0) {
  if (opt0->a) {  // matching score is changed
    if (!opt0->b) opt->b *= opt->a;
    if (!opt0->T) opt->T *= opt->a;
    if (!opt0->o_del) opt->o_del *= opt->a;
    if (!opt0->e_del) opt->e_del *= opt->a;
    if (!opt0->o_ins) opt->o_ins *= opt->a;
    if (!opt0->e_ins) opt->e_ins *= opt->a;
    if (!opt0->zdrop) opt->zdrop *= opt->a;
    if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
    if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
    if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
  }
}

// Manual option parsing to avoid getopt global state issues in a library
mem_opt_t* get_opts(int argc, char* argv[], bwaidx_t* idx) {
  mem_opt_t *opt, opt0;
  int i;
  char* p;

  bwa_verbose = 0;
  opt = mem_opt_init();
  memset(&opt0, 0, sizeof(mem_opt_t));

  for (i = 1; i < argc; ++i) {
    if (argv[i][0] != '-' || argv[i][1] == '\0') continue;
    char c = argv[i][1];
    
    // Simple manual parser for the subset of options used by the Python library
    if (c == 'k') { if (i+1 < argc) { opt->min_seed_len = atoi(argv[++i]); opt0.min_seed_len = 1; } }
    else if (c == 'c') { if (i+1 < argc) { opt->max_occ = atoi(argv[++i]); opt0.max_occ = 1; } }
    else if (c == 'T') { if (i+1 < argc) { opt->T = atoi(argv[++i]); opt0.T = 1; } }
    else if (c == 'U') { if (i+1 < argc) { opt->pen_unpaired = atoi(argv[++i]); opt0.pen_unpaired = 1; } }
    else if (c == 'a') opt->flag |= MEM_F_ALL;
    else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
    else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
    else if (c == 'v') { if (i+1 < argc) bwa_verbose = atoi(argv[++i]); }
    else if (c == 'L') { 
        if (i+1 < argc) {
            char *s = argv[++i];
            opt->pen_clip5 = opt->pen_clip3 = (int)strtol(s, &p, 10);
            if (*p != 0) opt->pen_clip3 = (int)strtol(p + 1, &p, 10);
            opt0.pen_clip5 = opt0.pen_clip3 = 1;
        }
    }
    else if (c == 'O') {
        if (i+1 < argc) {
            char *s = argv[++i];
            opt->o_del = opt->o_ins = (int)strtol(s, &p, 10);
            if (*p != 0) opt->o_ins = (int)strtol(p + 1, &p, 10);
            opt0.o_del = opt0.o_ins = 1;
        }
    }
    else if (c == 'E') {
        if (i+1 < argc) {
            char *s = argv[++i];
            opt->e_del = opt->e_ins = (int)strtol(s, &p, 10);
            if (*p != 0) opt->e_ins = (int)strtol(p + 1, &p, 10);
            opt0.e_del = opt0.e_ins = 1;
        }
    }
  }

  update_a(opt, &opt0);
  bwa_fill_scmat(opt->a, opt->b, opt->mat);

  if (idx && idx->bns) {
      // Handle any index-specific adjustments here if needed
  }

  return opt;
}
