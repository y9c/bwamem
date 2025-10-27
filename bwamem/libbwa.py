import argparse
from collections import namedtuple
import importlib
import importlib.util
import os
import sys  # noqa: F401
from contextlib import contextmanager

from cffi import FFI

ffi = FFI()

"""High-level interface to bwa (mem) aligner."""


@contextmanager
def suppress_stderr():
    """Context manager to suppress stderr output."""
    # Save the original stderr file descriptor
    old_stderr_fd = os.dup(2)
    try:
        # Redirect stderr to /dev/null
        with open(os.devnull, "w") as devnull:
            os.dup2(devnull.fileno(), 2)
            yield
    finally:
        # Restore the original stderr
        os.dup2(old_stderr_fd, 2)
        os.close(old_stderr_fd)


def get_shared_lib(name):
    """Cross-platform resolution of shared-object libraries, working
    around vagueries of setuptools.
    :param name: name of shared library to find.

    :returns: FFI shared library object.
    """
    try:
        # after 'python setup.py install' we should be able to do this
        lib_file = importlib.import_module(name).__file__
    except Exception:
        try:
            # after 'python setup.py develop' this should work
            spec = importlib.util.find_spec(name)
            if spec is None or spec.origin is None:
                raise ImportError(f'Cannot locate C library "{name}".')
            lib_file = spec.origin
        except Exception:
            raise ImportError('Cannot locate C library "{}".'.format(name))
        else:
            lib_file = os.path.abspath(lib_file)
    finally:
        library = ffi.dlopen(lib_file)
    return library


libbwa = get_shared_lib("bwalib")
# Reduce noisy internal BWA logs (mem_pestat notices) for tiny demos
try:
    libbwa.bwa_verbose = 0
except Exception:
    pass

ffi.cdef("""
  ////////////////////////////////
  // Alignment hit list structures
  //
  typedef struct {
    int64_t rb, re; // [rb,re): reference sequence in the alignment
    int qb, qe;     // [qb,qe): query sequence in the alignment
    int rid;        // reference seq ID
    int score;      // best local SW score
    int truesc;     // actual score corresponding to the aligned region; possibly smaller than $score
    int sub;        // 2nd best SW score
    int alt_sc;
    int csub;       // SW score of a tandem hit
    int sub_n;      // approximate number of suboptimal hits
    int w;          // actual band width used in extension
    int seedcov;    // length of regions coverged by seeds
    int secondary;  // index of the parent hit shadowing the current hit; <0 if primary
    int secondary_all;
    int seedlen0;   // length of the starting seed
    int n_comp:30, is_alt:2; // number of sub-alignments chained together
    float frac_rep;
    uint64_t hash;
  } mem_alnreg_t;

  typedef struct { size_t n, m; mem_alnreg_t *a; } mem_alnreg_v;

  typedef struct {   // This struct is only used for the convenience of API.
    int64_t pos;     // forward strand 5'-end mapping position
    int rid;         // reference sequence index in bntseq_t; <0 for unmapped
    int flag;        // extra flag
    uint32_t is_rev:1, is_alt:1, mapq:8, NM:22; // is_rev: whether on the reverse strand; mapq: mapping quality; NM: edit distance
    int n_cigar;     // number of CIGAR operations
    uint32_t *cigar; // CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234
    char *XA;        // alternative mappings

    int score, sub, alt_sc;
  } mem_aln_t;

  typedef struct { size_t n; mem_aln_t *aln; } mem_aln_v;

  void free_mem_aln_v (mem_aln_v *alns);

  ///////////////////////
  // Paired-end statistics
  //
  typedef struct {
    int failed;
    int low, high;
    double avg, std;
  } mem_pestat_t;

  ///////////////////////
  // Sequence structure
  //
  typedef struct {
    int l_seq, id;
    char *name, *comment, *seq, *qual, *sam;
  } bseq1_t;

  ///////////////////////
  // bwa index structures
  //
  typedef uint64_t bwtint_t;

  typedef struct {
    bwtint_t primary; // S^{-1}(0), or the primary index of BWT
    bwtint_t L2[5]; // C(), cumulative count
    bwtint_t seq_len; // sequence length
    bwtint_t bwt_size; // size of bwt, about seq_len/4
    uint32_t *bwt; // BWT
    // occurance array, separated to two parts
    uint32_t cnt_table[256];
    // suffix array
    int sa_intv;
    bwtint_t n_sa;
    bwtint_t *sa;
  } bwt_t;

  typedef struct {
    int64_t offset;
    int32_t len;
    int32_t n_ambs;
    uint32_t gi;
    int32_t is_alt;
    char *name, *anno;
  } bntann1_t;

  typedef struct {
    int64_t offset;
    int32_t len;
    char amb;
  } bntamb1_t;

  typedef struct {
    int64_t l_pac;
    int32_t n_seqs;
    uint32_t seed;
    bntann1_t *anns; // n_seqs elements
    int32_t n_holes;
    bntamb1_t *ambs; // n_holes elements
    FILE *fp_pac;
  } bntseq_t;

  typedef struct {
    bwt_t    *bwt; // FM-index
    bntseq_t *bns; // information on the reference sequences
    uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base

    int    is_shm;
    int64_t l_mem;
    uint8_t  *mem;
  } bwaidx_t;

  bwaidx_t *bwa_idx_load_all(const char *hint);
  void bwa_idx_destroy(bwaidx_t *idx);

  /////////////////
  // Option parsing
  //
  typedef struct {
    int a, b;               // match score and mismatch penalty
    int o_del, e_del;
    int o_ins, e_ins;
    int pen_unpaired;       // phred-scaled penalty for unpaired reads
    int pen_clip5,pen_clip3;// clipping penalty. This score is not deducted from the DP score.
    int w;                  // band width
    int zdrop;              // Z-dropoff

    uint64_t max_mem_intv;

    int T;                  // output score threshold; only affecting output
    int flag;               // see MEM_F_* macros
    int min_seed_len;       // minimum seed length
    int min_chain_weight;
    int max_chain_extend;
    float split_factor;     // split into a seed if MEM is longer than min_seed_len*split_factor
    int split_width;        // split into a seed if its occurence is smaller than this value
    int max_occ;            // skip a seed if its occurence is larger than this value
    int max_chain_gap;      // do not chain seed if it is max_chain_gap-bp away from the closest seed
    int n_threads;          // number of threads
    int chunk_size;         // process chunk_size-bp sequences in a batch
    float mask_level;       // regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
    float drop_ratio;       // drop a chain if its seed coverage is below drop_ratio times the seed coverage of a better chain overlapping with the small chain
    float XA_drop_ratio;    // when counting hits for the XA tag, ignore alignments with score < XA_drop_ratio * max_score; only effective for the XA tag
    float mask_level_redun;
    float mapQ_coef_len;
    int mapQ_coef_fac;
    int max_ins;            // when estimating insert size distribution, skip pairs with insert longer than this value
    int max_matesw;         // perform maximally max_matesw rounds of mate-SW for each end
    int max_XA_hits, max_XA_hits_alt; // if there are max_hits or fewer, output them all
    int8_t mat[25];         // scoring matrix; mat[0] == 0 if unset
  } mem_opt_t;

  mem_opt_t * get_opts(int argc, char *argv[], bwaidx_t * idx);

  static const char valid_opts[];


  ///////////////////
  // Run an alignment
  //
  mem_aln_v *align(mem_opt_t * opt, bwaidx_t * index, char * seq);

  ///////////////////
  // Single-end alignment functions
  //
  mem_alnreg_v mem_align1(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq);
  mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq, const mem_alnreg_t *ar);
  // Wrapper that returns a pointer to avoid bitfield-return limitation
  mem_aln_t *mem_reg2aln_ptr(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq, const mem_alnreg_t *ar);

  ///////////////////
  // Paired-end alignment functions
  //
  void mem_pestat(const mem_opt_t *opt, int64_t l_pac, int n, const mem_alnreg_v *regs, mem_pestat_t pes[4]);
  int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2]);

  ///////////////////
  // Index building
  //
  int bwa_idx_build(const char *fa, const char *prefix, int algo_type, int block_size);
  
  ///////////////////
  // FASTQ reading functions
  //
  // kstring_t structure
  typedef struct {
    size_t l, m;
    char *s;
  } kstring_t;
  
  // kseq_t structure
  typedef struct {
    kstring_t name, comment, seq, qual;
    int last_char;
    void *f;
  } kseq_t;
  
  // File opening functions
  void *err_xzopen_core(const char *func, const char *fn, const char *mode);
  int gzclose(void *file);
  // kseq functions for FASTQ reading
  kseq_t *kseq_init(void *fp);
  int kseq_read(kseq_t *ks);
  void kseq_destroy(kseq_t *ks);
  // bseq functions for batch reading
  bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_);
  
  // C standard free for releasing memory allocated by BWA
  void free(void *ptr);
  // Control BWA logging verbosity
  extern int bwa_verbose;
""")


# Alignment result for single reads
Alignment = namedtuple(
    "Alignment",
    ["rname", "orient", "pos", "mapq", "cigar", "NM", "score", "is_primary"],
)

# Paired-end alignment result
PairedAlignment = namedtuple(
    "PairedAlignment", ["read1", "read2", "is_proper_pair", "insert_size"]
)


class BwaAligner(object):
    def __init__(
        self,
        index: str,
        options: str = "",
        *,
        min_score: int | None = None,
        softclip_supplementary: bool | None = None,
        mark_secondary: bool | None = None,
        clip_penalties: tuple[int, int] | None = None,
        unpaired_penalty: int | None = None,
        insert_model: tuple | None = None,
    ):
        """Interface to bwa mem alignment.

        :param index: bwa index base path.
        :param options: alignment options as would be given on the bwa mem command line.
        :param min_score: equivalent to -T (minimum alignment score to output)
        :param softclip_supplementary: equivalent to -Y (soft-clip supplementary alignments)
        :param mark_secondary: equivalent to -M (mark shorter split hits as secondary)
        :param clip_penalties: equivalent to -L a,b (5'/3' clipping penalties)
        :param unpaired_penalty: equivalent to -U (penalty for unpaired alignments)
        :param insert_model: equivalent to -I mean[,std[,max[,min]]] for FR pairing model

        """
        self.index_base = index.encode()
        self._cigchar = "MIDSH"

        # Compose CLI-style options from explicit kwargs
        extra_opts: list[str] = []
        if softclip_supplementary:
            extra_opts += ["-Y"]
        if mark_secondary:
            extra_opts += ["-M"]
        if clip_penalties is not None:
            a, b = clip_penalties
            extra_opts += ["-L", f"{int(a)},{int(b)}"]
        if unpaired_penalty is not None:
            extra_opts += ["-U", str(int(unpaired_penalty))]
        if min_score is not None:
            extra_opts += ["-T", str(int(min_score))]
        # Support -I mean,std,max[,min]
        self._insert_model = None
        if insert_model is not None:
            n = len(insert_model)
            if n < 1 or n > 4:
                raise ValueError(
                    "insert_model must be (mean), (mean,std), (mean,std,max), or (mean,std,max,min)"
                )
            mean = float(insert_model[0])
            std = float(insert_model[1]) if n >= 2 else max(1.0, mean * 0.1)
            imax = int(insert_model[2]) if n >= 3 else int(mean + 4.0 * std + 0.499)
            imin = int(insert_model[3]) if n >= 4 else None
            # Build -I string respecting provided arity
            if n == 1:
                i_arg = f"{mean}"
            elif n == 2:
                i_arg = f"{mean},{std}"
            elif n == 3:
                i_arg = f"{mean},{std},{imax}"
            else:
                i_arg = f"{mean},{std},{imax},{imin}"
            extra_opts += ["-I", i_arg]
            # Store normalized model for pairing logic
            self._insert_model = (mean, std, imax) + (
                (imin,) if imin is not None else tuple()
            )

        # Merge explicit kwargs options with provided options string
        merged_options = " ".join(extra_opts + (options.split() if options else []))

        # Validate options against allowed set
        valid_opts = ffi.string(libbwa.valid_opts).decode().replace(":", "")
        for opt in merged_options.split():
            if opt[0] != "-":
                continue
            if opt[1] not in valid_opts:
                raise ValueError(
                    "Option '{}' is not a valid option (allowed: {}).".format(
                        opt, " ".join(valid_opts)
                    )
                )

        # we need to pass the index to the option parsing
        # TODO: clean up this requirement
        self.index = libbwa.bwa_idx_load_all(self.index_base)
        if self.index == ffi.NULL:
            raise ValueError("Failed to load bwa index.")

        argv = ["bwamem"] + (merged_options.split() if merged_options else [])
        argc = len(argv)
        self.opt = libbwa.get_opts(
            argc, [ffi.new("char[]", x.encode()) for x in argv], self.index
        )
        if self.opt == ffi.NULL:
            raise ValueError("Failed to parse options.")

    def __del__(self):
        if hasattr(self, "index"):
            try:
                libbwa.bwa_idx_destroy(self.index)
            except (AttributeError, NameError):
                # Function not available, skip cleanup
                pass
        # No additional opt memory to free

    def align(self, seq1: str, seq2: str = None):
        """Align one or two sequences to the index.

        :param seq1: first sequence to align (required)
        :param seq2: second sequence to align (optional, for paired-end)

        :returns: If seq2 is None, returns tuple of :class:`Alignment` for single-end.
                 If seq2 is provided, returns tuple of :class:`PairedAlignment` for paired-end.
        """
        if seq2 is None:
            # Single-end alignment
            return self._align_single_end(seq1)
        else:
            # Paired-end alignment
            return self._align_paired_end(seq1, seq2)

    def _align_single_end(self, seq: str):
        """Perform single-end alignment using the new BWA functions."""
        # Get alignment regions
        regs = libbwa.mem_align1(
            self.opt,
            self.index.bwt,
            self.index.bns,
            self.index.pac,
            len(seq),
            seq.encode(),
        )

        if regs.n == 0:
            return tuple()

        # Convert regions to alignments
        alignments = []
        for i in range(regs.n):
            if regs.a[i].score >= self.opt.T:  # Only keep alignments above threshold
                aln_ptr = libbwa.mem_reg2aln_ptr(
                    self.opt,
                    self.index.bns,
                    self.index.pac,
                    len(seq),
                    seq.encode(),
                    ffi.addressof(regs.a[i]),
                )
                if aln_ptr != ffi.NULL and aln_ptr.rid >= 0:  # Valid alignment
                    cigar = self._build_cigar(aln_ptr.cigar, aln_ptr.n_cigar)
                    alignment = Alignment(
                        rname=ffi.string(
                            self.index.bns.anns[aln_ptr.rid].name
                        ).decode(),
                        orient="+-"[aln_ptr.is_rev],
                        pos=aln_ptr.pos,
                        mapq=aln_ptr.mapq,
                        cigar=cigar,
                        NM=aln_ptr.NM,
                        score=regs.a[i].score,
                        is_primary=(i == 0),
                    )
                    alignments.append(alignment)
                    # Free dynamically allocated CIGAR and struct
                    if aln_ptr.cigar != ffi.NULL:
                        libbwa.free(aln_ptr.cigar)
                    libbwa.free(aln_ptr)

        # Free the mem_alnreg_v array allocated by mem_align1
        if regs.a != ffi.NULL:
            libbwa.free(regs.a)

        return tuple(alignments)

    def _align_paired_end(self, seq1: str, seq2: str):
        """Perform paired-end alignment."""
        # If an insert_model was set on the aligner, use it to override per-call args
        eff_insert_size = None
        eff_insert_std = None
        eff_insert_min = None
        eff_insert_max = None
        if getattr(self, "_insert_model", None) is not None:
            model = self._insert_model
            # model = (mean, std, max) or (mean, std, max, min)
            eff_insert_size = float(model[0])
            eff_insert_std = float(model[1])
            eff_insert_max = int(model[2])
            if len(model) >= 4:
                eff_insert_min = int(model[3])
        # Get alignment regions for both reads
        regs1 = libbwa.mem_align1(
            self.opt,
            self.index.bwt,
            self.index.bns,
            self.index.pac,
            len(seq1),
            seq1.encode(),
        )
        regs2 = libbwa.mem_align1(
            self.opt,
            self.index.bwt,
            self.index.bns,
            self.index.pac,
            len(seq2),
            seq2.encode(),
        )

        # Create arrays for paired-end processing
        regs_array = ffi.new("mem_alnreg_v[2]")
        regs_array[0] = regs1
        regs_array[1] = regs2

        # Set up insert size distribution
        pes = ffi.new("mem_pestat_t[4]")
        if eff_insert_size is not None:
            # Use provided insert size
            pes[1].failed = 0
            pes[1].avg = eff_insert_size
            pes[1].std = (
                eff_insert_std if eff_insert_std is not None else eff_insert_size * 0.1
            )
            pes[1].high = int(pes[1].avg + 4.0 * pes[1].std + 0.499)
            pes[1].low = int(pes[1].avg - 4.0 * pes[1].std + 0.499)
            if pes[1].low < 1:
                pes[1].low = 1
            # Apply explicit min/max overrides if given
            if eff_insert_max is not None:
                pes[1].high = int(eff_insert_max)
            if eff_insert_min is not None:
                pes[1].low = int(eff_insert_min)
        else:
            # Infer insert size from data
            with suppress_stderr():
                libbwa.mem_pestat(self.opt, self.index.bns.l_pac, 2, regs_array, pes)

        # Skip mem_sam_pe; construct paired alignments from regions

        # Parse the results (this is simplified - in practice you'd parse the SAM output)
        # For now, return basic paired-end information
        paired_alignments = []

        # Convert regions to single alignments for each read
        read1_alignments = self._convert_regions_to_alignments(regs1, seq1, 1)
        read2_alignments = self._convert_regions_to_alignments(regs2, seq2, 2)
        # Free mem_alnreg arrays
        if regs1.a != ffi.NULL:
            libbwa.free(regs1.a)
        if regs2.a != ffi.NULL:
            libbwa.free(regs2.a)

        # Create paired alignments
        for aln1 in read1_alignments:
            for aln2 in read2_alignments:
                # Check if this is a proper pair
                is_proper = self._is_proper_pair(
                    aln1,
                    aln2,
                    pes,
                    len(seq1),
                    len(seq2),
                    eff_insert_size,
                    eff_insert_std,
                    eff_insert_min,
                    eff_insert_max,
                )
                insert_size_val = (
                    self._calculate_insert_size(aln1, aln2, len(seq1), len(seq2))
                    if is_proper
                    else None
                )

                paired_aln = PairedAlignment(
                    read1=aln1,
                    read2=aln2,
                    is_proper_pair=is_proper,
                    insert_size=insert_size_val,
                )
                paired_alignments.append(paired_aln)

        return tuple(paired_alignments)

    def _convert_regions_to_alignments(self, regs, seq, read_num):
        """Convert alignment regions to Alignment objects."""
        alignments = []
        for i in range(regs.n):
            if regs.a[i].score >= self.opt.T:
                aln_ptr = libbwa.mem_reg2aln_ptr(
                    self.opt,
                    self.index.bns,
                    self.index.pac,
                    len(seq),
                    seq.encode(),
                    ffi.addressof(regs.a[i]),
                )
                if aln_ptr != ffi.NULL and aln_ptr.rid >= 0:
                    cigar = self._build_cigar(aln_ptr.cigar, aln_ptr.n_cigar)
                    alignment = Alignment(
                        rname=ffi.string(
                            self.index.bns.anns[aln_ptr.rid].name
                        ).decode(),
                        orient="+-"[aln_ptr.is_rev],
                        pos=aln_ptr.pos,
                        mapq=aln_ptr.mapq,
                        cigar=cigar,
                        NM=aln_ptr.NM,
                        score=regs.a[i].score,
                        is_primary=(i == 0),
                    )
                    alignments.append(alignment)
                    if aln_ptr.cigar != ffi.NULL:
                        libbwa.free(aln_ptr.cigar)
                    libbwa.free(aln_ptr)
        return alignments

    def _is_proper_pair(
        self,
        aln1,
        aln2,
        pes,
        len1: int,
        len2: int,
        user_insert: float | None,
        user_std: float | None,
        user_min: int | None,
        user_max: int | None,
    ):
        """Check if two alignments form a proper FR pair with plausible insert size."""
        if aln1.rname != aln2.rname:
            return False

        # Require FR orientation: read1 on '+' before read2 on '-'
        is_fr = (aln1.orient == "+") and (aln2.orient == "-") and (aln1.pos <= aln2.pos)
        if not is_fr:
            return False

        insert_size_val = self._calculate_insert_size(aln1, aln2, len1, len2)
        if insert_size_val is None:
            return False

        # Build candidate window
        low = None
        high = None

        # Explicit min/max override everything if provided
        if user_min is not None:
            low = int(user_min)
        if user_max is not None:
            high = int(user_max)

        # If explicit bounds are incomplete, fill from avg/std if available
        if (low is None or high is None) and user_insert is not None:
            std = user_std if user_std is not None else max(1.0, user_insert * 0.1)
            avg_low = int(user_insert - 4.0 * std + 0.499)
            avg_high = int(user_insert + 4.0 * std + 0.499)
            if avg_low < 1:
                avg_low = 1
            if low is None:
                low = avg_low
            if high is None:
                high = avg_high

        # If still missing, attempt to use FR bin from pes
        if low is None or high is None:
            try:
                fr_pes = pes[1]
                if not fr_pes.failed:
                    if low is None:
                        low = int(fr_pes.low)
                    if high is None:
                        high = int(fr_pes.high)
            except Exception:
                pass

        # If we obtained any bounds, enforce them strictly
        if low is not None or high is not None:
            if low is None:
                return insert_size_val <= high
            if high is None:
                return insert_size_val >= low
            return low <= insert_size_val <= high

        # No bounds available: reject to be strict
        return False

    def _calculate_insert_size(self, aln1, aln2, len1: int, len2: int):
        """Calculate insert size for FR pairs (distance between 5' ends including read2 length)."""
        if aln1.rname != aln2.rname:
            return None
        if aln1.orient == "+" and aln2.orient == "-" and aln1.pos <= aln2.pos:
            # 5' ends: aln1.pos and aln2.pos + len2 - 1; include read2 length
            return int((aln2.pos + len2) - aln1.pos)
        return None

    def _build_cigar(self, cigar_array, n_cigar):
        """Build CIGAR string from CIGAR array."""
        if n_cigar == 0:
            return ""

        cigar = ""
        for i in range(n_cigar):
            op_len = cigar_array[i] >> 4
            op = cigar_array[i] & 0xF
            cigar += str(op_len) + self._cigchar[op]

        return cigar


class BwaIndexer(object):
    """Interface to BWA index building functionality."""

    # Algorithm types
    BWTALGO_AUTO = 0
    BWTALGO_RB2 = 1
    BWTALGO_BWTSW = 2
    BWTALGO_IS = 3

    def __init__(self, algorithm="auto", block_size=10000000):
        """Initialize BWA indexer.

        :param algorithm: BWT construction algorithm ('auto', 'rb2', 'bwtsw', 'is')
        :param block_size: Block size for bwtsw algorithm (in bytes)
        """
        self.algorithm = algorithm
        self.block_size = block_size

        # Convert algorithm string to integer
        algo_map = {
            "auto": self.BWTALGO_AUTO,
            "rb2": self.BWTALGO_RB2,
            "bwtsw": self.BWTALGO_BWTSW,
            "is": self.BWTALGO_IS,
        }

        if algorithm not in algo_map:
            raise ValueError(
                f"Unknown algorithm '{algorithm}'. Choose from: {list(algo_map.keys())}"
            )

        self.algo_type = algo_map[algorithm]

    def build_index(self, fasta_file, prefix=None):
        """Build BWA index from FASTA file.

        :param fasta_file: Path to input FASTA file
        :param prefix: Output prefix for index files (default: same as FASTA file)
        :returns: Path to the index prefix
        """
        import os

        if not os.path.exists(fasta_file):
            raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

        if prefix is None:
            # Use FASTA filename without extension as prefix
            prefix = os.path.splitext(fasta_file)[0]

        # Convert to bytes for C function
        fasta_bytes = fasta_file.encode("utf-8")
        prefix_bytes = prefix.encode("utf-8")

        # Call the C function directly
        result = libbwa.bwa_idx_build(
            fasta_bytes, prefix_bytes, self.algo_type, self.block_size
        )

        if result != 0:
            raise RuntimeError(f"Failed to build BWA index for {fasta_file}")

        return prefix

    def build_index_with_options(
        self, fasta_file, prefix=None, algorithm=None, block_size=None
    ):
        """Build BWA index with specific options.

        :param fasta_file: Path to input FASTA file
        :param prefix: Output prefix for index files
        :param algorithm: BWT construction algorithm ('auto', 'rb2', 'bwtsw', 'is')
        :param block_size: Block size for bwtsw algorithm
        :returns: Path to the index prefix
        """
        # Temporarily override instance settings
        old_algorithm = self.algorithm
        old_block_size = self.block_size

        if algorithm is not None:
            self.algorithm = algorithm
            algo_map = {
                "auto": self.BWTALGO_AUTO,
                "rb2": self.BWTALGO_RB2,
                "bwtsw": self.BWTALGO_BWTSW,
                "is": self.BWTALGO_IS,
            }
            if algorithm not in algo_map:
                raise ValueError(
                    f"Unknown algorithm '{algorithm}'. "
                    f"Choose from: {list(algo_map.keys())}"
                )
            self.algo_type = algo_map[algorithm]

        if block_size is not None:
            self.block_size = block_size

        try:
            result = self.build_index(fasta_file, prefix)
        finally:
            # Restore original settings
            self.algorithm = old_algorithm
            self.block_size = old_block_size
            if algorithm is not None:
                algo_map = {
                    "auto": self.BWTALGO_AUTO,
                    "rb2": self.BWTALGO_RB2,
                    "bwtsw": self.BWTALGO_BWTSW,
                    "is": self.BWTALGO_IS,
                }
                self.algo_type = algo_map[self.algorithm]

        return result


def get_parser():
    parser = argparse.ArgumentParser("Align a sequence with bwa mem.")
    parser.add_argument("index", help="bwa index base path.")
    parser.add_argument("sequence", nargs="+", help="base sequence")
    return parser


def main():
    args, opts = get_parser().parse_known_args()
    options = ""
    if len(opts) > 0:
        options = " ".join(opts)
    aligner = BwaAligner(args.index, options=options)
    for i, seq in enumerate(args.sequence, 1):
        alignments = aligner.align_seq(seq)
        print("Found {} alignments for input {}.".format(len(alignments), i))
        for aln in alignments:
            print("  ", aln)
