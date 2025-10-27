#!/usr/bin/env python3
"""Create demo FASTA and FASTQ files for testing."""

def create_demo_files():
    """Create demo FASTA and FASTQ files for testing (under tests/)."""
    import os

    os.makedirs('tests', exist_ok=True)

    # Create demo FASTA reference
    fasta_content = ">chr1\n" + ("ATCG" * 64) + "\n" + ("ATCG" * 64) + "\n"
    with open('tests/demo_reference.fasta', 'w') as f:
        f.write(fasta_content)

    # Create demo FASTQ R1 file with matching seq/qual lengths (64)
    r1_1 = "A" * 64
    r1_2 = "C" * 64
    r1_3 = "G" * 64
    q1_1 = "I" * 64
    q1_2 = "H" * 64
    q1_3 = "G" * 64
    fastq_r1_content = f"""@read1/1
{r1_1}
+
{q1_1}
@read2/1
{r1_2}
+
{q1_2}
@read3/1
{r1_3}
+
{q1_3}
"""
    with open('tests/demo_reads_R1.fastq', 'w') as f:
        f.write(fastq_r1_content)

    # Create demo FASTQ R2 file with matching seq/qual lengths (64)
    r2_1 = "T" * 64
    r2_2 = "G" * 64
    r2_3 = "C" * 64
    q2_1 = "I" * 64
    q2_2 = "H" * 64
    q2_3 = "G" * 64
    fastq_r2_content = f"""@read1/2
{r2_1}
+
{q2_1}
@read2/2
{r2_2}
+
{q2_2}
@read3/2
{r2_3}
+
{q2_3}
"""
    with open('tests/demo_reads_R2.fastq', 'w') as f:
        f.write(fastq_r2_content)

    print("Created demo files:")
    print("- tests/demo_reference.fasta")
    print("- tests/demo_reads_R1.fastq")
    print("- tests/demo_reads_R2.fastq")

if __name__ == "__main__":
    create_demo_files()
