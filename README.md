bwamem
======

Python bindings to `bwa mem` aligner; sufficient to load and index and perform
alignments of sequences to the index to obtain basic statistics.

These python bindings are licensed under Mozilla Public License 2.0, bwa is licenced
under GNU General Public License v3.0.

Documentation can be found at https://y9c.github.io/bwamem/.

Installation
------------

The git source repository contains bwa as a submodule. The repository should therefore
be cloned using the recursive option.

The package `setup.py` script requires `libbwa.a` to have been built in the submodule
directory before running. This can be performed via the `libbwa.a` target, which first
makes some amendments to the bwa/Makefile. To build and install the package one should
therefore run:

    git clone --recursive https://github.com/y9c/bwamem.git
    cd bwamem
    make bwa/libbwa.a 
    python setup.py install


Building BWA Indexes
--------------------

The `BwaIndexer` class provides a pythonic interface to build BWA indexes from
FASTA files. It supports different BWT construction algorithms:

```python
from bwamem import BwaIndexer

# Create indexer with default settings (auto algorithm)
indexer = BwaIndexer()

# Build index from FASTA file
index_path = indexer.build_index('reference.fa')
print(f'Index built at: {index_path}')

# Use specific algorithm
indexer = BwaIndexer(algorithm='is')  # or 'rb2', 'bwtsw', 'auto'
index_path = indexer.build_index('reference.fa', prefix='my_index')
```

Available algorithms:
- `auto`: Automatically choose algorithm based on genome size
- `rb2`: RB2 algorithm (good for medium genomes)
- `bwtsw`: BWT-SW algorithm (good for large genomes)
- `is`: IS algorithm (good for small genomes)

Performing Alignments
---------------------

The `BwaAligner` class provides a pythonic interface to `bwa mem` aligner. It
takes as input a bwa index fileset on construction and can then be used to find
alignments of sequences given as strings.

### Single-End Alignment

For single-end reads, use the `align()` method with one sequence:

```python
from bwamem import BwaAligner, Alignment
index = 'path/to/index' # the path given to bwa index
seq = 'ACGATCGCGATCGA'

aligner = BwaAligner(index)
alignments = aligner.align(seq)  # Returns tuple of Alignment objects
print('Found {} alignments.'.format(len(alignments)))
for aln in alignments:
    print(f'  {aln.rname}:{aln.pos} {aln.orient} {aln.cigar} (mapq={aln.mapq}, score={aln.score})')
```

### Paired-End Alignment

For paired-end reads, use the `align()` method with two sequences:

```python
from bwamem import BwaAligner, PairedAlignment
index = 'path/to/index'
read1 = 'ACGATCGCGATCGA'
read2 = 'TTCGATCGATCGAT'

aligner = BwaAligner(index)
paired_alignments = aligner.align(read1, read2)  # Returns tuple of PairedAlignment objects
print('Found {} paired alignments.'.format(len(paired_alignments)))
for pe_aln in paired_alignments:
    print(f'  Read1: {pe_aln.read1.rname}:{pe_aln.read1.pos} {pe_aln.read1.orient}')
    print(f'  Read2: {pe_aln.read2.rname}:{pe_aln.read2.pos} {pe_aln.read2.orient}')
    print(f'  Proper pair: {pe_aln.is_proper_pair}, Insert size: {pe_aln.insert_size}')
```

### Custom Insert Size Distribution

For paired-end reads, you can specify the expected insert size distribution:

```python
# With custom insert size parameters
paired_alignments = aligner.align(read1, read2, insert_size=500, insert_std=50)
```

### Data Structures

**Alignment** (for single-end reads):
```python
Alignment(rname='chr1', orient='+', pos=1000, mapq=60, cigar='100M', NM=0, score=100, is_primary=True)
```

**PairedAlignment** (for paired-end reads):
```python
PairedAlignment(read1=Alignment(...), read2=Alignment(...), is_proper_pair=True, insert_size=500)
```

### Alignment Parameters

Alignment parameters can be given as they are on the `bwa mem` command line:

```python
from bwamem import BwaAligner
index = 'path/to/index'
options = '-x ont2d -A 1 -B 0'
aligner = BwaAligner(index, options=options)
```

The package now supports all BWA MEM options including paired-end specific parameters like insert size distribution (`-I` option).

Complete Workflow Example
-------------------------

Here's a complete example showing how to build an index and perform both single-end and paired-end alignments:

```python
from bwamem import BwaIndexer, BwaAligner, Alignment, PairedAlignment

# Step 1: Build index from FASTA file
indexer = BwaIndexer(algorithm='auto')
index_path = indexer.build_index('reference.fa')
print(f'Index built at: {index_path}')

# Step 2: Create aligner with the index
aligner = BwaAligner(index_path)

# Step 3a: Single-end alignment
print("Single-end alignments:")
se_sequences = ['ACGATCGCGATCGA', 'GCTAGCTAGCTAG']
for i, seq in enumerate(se_sequences, 1):
    alignments = aligner.align(seq)
    print(f'Found {len(alignments)} alignments for sequence {i}')
    for aln in alignments:
        print(f'  {aln.rname}:{aln.pos} {aln.orient} {aln.cigar} (mapq={aln.mapq})')

# Step 3b: Paired-end alignment
print("\nPaired-end alignments:")
pe_reads = [
    ('ACGATCGCGATCGA', 'TTCGATCGATCGAT'),
    ('GCTAGCTAGCTAG', 'CGATCGATCGATC')
]
for i, (read1, read2) in enumerate(pe_reads, 1):
    paired_alignments = aligner.align(read1, read2)
    print(f'Found {len(paired_alignments)} paired alignments for read pair {i}')
    for pe_aln in paired_alignments:
        print(f'  Read1: {pe_aln.read1.rname}:{pe_aln.read1.pos} {pe_aln.read1.orient}')
        print(f'  Read2: {pe_aln.read2.rname}:{pe_aln.read2.pos} {pe_aln.read2.orient}')
        print(f'  Proper pair: {pe_aln.is_proper_pair}, Insert size: {pe_aln.insert_size}')

# Step 3c: Paired-end with custom insert size
print("\nPaired-end with custom insert size:")
pe_alignments = aligner.align('ACGATCGCGATCGA', 'TTCGATCGATCGAT', 
                             insert_size=500, insert_std=50)
print(f'Found {len(pe_alignments)} alignments with custom insert size')
```

### Advanced Paired-End Features

```python
# Filter for proper pairs only
proper_pairs = [pe for pe in pe_alignments if pe.is_proper_pair]
print(f'Found {len(proper_pairs)} proper pairs')

# Access individual read alignments
for pe_aln in pe_alignments:
    read1_aln = pe_aln.read1
    read2_aln = pe_aln.read2
    if read1_aln.is_primary and read2_aln.is_primary:
        print(f'Primary alignment: {read1_aln.rname}:{read1_aln.pos}-{read2_aln.pos}')

```
