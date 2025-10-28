bwamem
======

Python bindings for the BWA-MEM aligner.

Installation
------------

```bash
pip install bwamem
```

Usage
-----

### Build Index

```python
from bwamem import BwaIndexer

indexer = BwaIndexer()
index_path = indexer.build_index('reference.fa')
```

### Single-End Alignment

```python
from bwamem import BwaAligner

aligner = BwaAligner('path/to/index')
alignments = aligner.align('ACGATCGCGATCGA')

for aln in alignments:
    print(f'{aln.ctg}:{aln.r_st} strand={aln.strand} mapq={aln.mapq}')
```

### Paired-End Alignment

```python
read1 = 'ACGATCGCGATCGA'
read2 = 'TTCGATCGATCGAT'

paired_alignments = aligner.align(read1, read2)

for pe_aln in paired_alignments:
    print(f'Insert size: {pe_aln.insert_size}, Proper pair: {pe_aln.is_proper_pair}')
```

### Custom Options

```python
# Specify alignment parameters
aligner = BwaAligner('path/to/index', options='-x ont2d -A 1 -B 0')

# Set custom insert size for paired-end reads
paired_alignments = aligner.align(read1, read2, insert_size=500, insert_std=50)
```

Alignment Attributes
--------------------

Each `Alignment` object contains the following attributes:

| Attribute | Description |
|-----------|-------------|
| `ctg` | Contig/reference name |
| `r_st` | Reference start position (0-based) |
| `r_en` | Reference end position (property) |
| `strand` | Strand: +1 for forward, -1 for reverse |
| `q_st`, `q_en` | Query start/end positions |
| `mapq` | Mapping quality score |
| `cigar` | CIGAR as list of `[length, op]` pairs |
| `cigar_str` | CIGAR string (property) |
| `NM` | Edit distance |
| `score` | Alignment score |
| `is_primary` | Primary alignment flag |

**Calculated properties** (computed on demand): `r_en`, `cigar_str`, `blen`, `mlen`

**CIGAR operations**: 0=M (match), 1=I (insertion), 2=D (deletion), 3=N (skip), 4=S (soft-clip), 5=H (hard-clip)

**PairedAlignment** contains: `read1`, `read2` (Alignment objects), `is_proper_pair` (bool), `insert_size` (int or None)

License
-------

- Python bindings: Mozilla Public License 2.0  
- BWA: GNU General Public License v3.0
