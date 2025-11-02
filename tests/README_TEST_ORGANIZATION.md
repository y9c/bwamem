# Test Organization

This document describes the organization of tests in the `tests/` directory.

## Test Files

### Automated pytest Tests

1. **`test_matesw.py`** - Tests for mem_matesw functionality and paired-end alignment
   - Tests alignment production
   - Tests rRNA hit detection (mate rescue)
   - Tests alignment structure and properties
   - Tests primary alignment marking
   - Tests same-contig pairing
   - Tests alignment scores

2. **`test_visualization.py`** - Tests for alignment visualization functionality
   - Tests `Alignment.visualize()` method
   - Tests visualization formatting and alignment
   - Tests `visualize_paired_alignment()` function
   - Tests pairing information display
   - Tests reverse complement handling
   - Tests line width parameter

3. **Other test files:**
   - `test_bwamem.py` - Basic BWA functionality tests
   - `test_cigar_insertion.py` - CIGAR operation tests
   - `test_fastx_reader.py` - FASTX file reading tests
   - `test_input_validation.py` - Input validation tests
   - `test_paired_mapping.py` - Paired-end mapping workflow tests

### Manual/Debug Scripts

- **`test_matesw_direct.py`** - Standalone script for manual testing/debugging
  - Provides detailed output and visualization
  - Useful for interactive debugging
  - Not run as part of pytest suite

## Running Tests

By default, pytest is configured to show all output (`-s` flag), so you'll see `print()` statements and test progress in real-time.

### Run all tests:
```bash
pytest tests/
```

### Run specific test files:
```bash
pytest tests/test_matesw.py
pytest tests/test_visualization.py
```

### Run with verbose output (default, already enabled):
```bash
pytest tests/test_matesw.py -v
pytest tests/test_visualization.py -v
```

### Run specific test functions:
```bash
pytest tests/test_matesw.py::test_matesw_rRNA_hits -v
pytest tests/test_visualization.py::test_visualize_paired_alignment_function -v
```

## Test Fixtures

### Shared Fixtures (in test files)

- `test_index_path` - Path to test reference index (`tests/test_data/reference/ref.mk.subset`)
- `aligner` - BwaAligner instance initialized with test index
- `paired_alignments` - Results from paired-end alignment of test sequences
- `rRNA_pair` - A paired alignment with rRNA hits (for visualization tests)

## Test Data

Test reference index files are located in:
- `tests/test_data/reference/ref.mk.subset.*` - Subset reference containing 79 sequences used in tests

Test sequences (MK converted):
- `SEQ1_CONV = "TTTTGGTTTTGGGTGGGGGTTGTTGGGGGGGGTGTTGTGGGGGTGGTT"`
- `SEQ2_CONV = "GGTTTGTGGTGGTTGGGTGTTTGTGGTGGTGTTGTTTTTTGGTTTTTTG"`

## Test Coverage

The tests cover:
- ✅ mem_matesw functionality
- ✅ Mate rescue (rRNA hit detection)
- ✅ Alignment structure and properties
- ✅ Visualization formatting
- ✅ Pairing information display
- ✅ Reverse complement handling
- ✅ Line width parameter handling

