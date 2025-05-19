# DNA Sequencing Algorithms Toolkit

A Python package for DNA sequence analysis, including FASTQ parsing, motif searching, quality score analysis, and GC content visualization.

## Features

- **FASTQ Parsing**: Parse FASTQ files to extract sequences and quality scores
- **Pattern Matching**: Find exact and approximate matches (with mismatches) of patterns within sequences
- **Reverse Complement**: Calculate reverse complements of DNA sequences
- **Quality Score Analysis**: Convert between ASCII and numeric quality scores, generate quality histograms
- **GC Content Analysis**: Calculate and visualize GC content distribution across sequence positions

## Installation

```bash
pip install git+https://github.com/andrewbrowne3/Algorithms_for_DNA_sequencing.git
```

## Usage Examples

### Reading FASTQ Files

```python
from dna_sequencing import fastq_parser

# Read sequences and quality scores from a FASTQ file
sequences, quality_scores = fastq_parser.readFastQ("sample.fastq")
print(f"Read {len(sequences)} sequences")
```

### Pattern Matching

```python
from dna_sequencing import fastq_parser

# Find all exact occurrences of a pattern in a text
pattern = "AGTC"
text = "AGTCAGTCAGTC"
matches = fastq_parser.naive(pattern, text)
print(f"Pattern found at positions: {matches}")

# Find all approximate matches with up to 2 mismatches
approx_matches = fastq_parser.naive_with_mismatches(pattern, text, max_mismatches=2)
print(f"Pattern found with mismatches at positions: {approx_matches}")
```

### Reverse Complement

```python
from dna_sequencing import fastq_parser

sequence = "AGTCN"
rev_comp = fastq_parser.ReverseComplement(sequence)
print(f"Reverse complement of {sequence} is {rev_comp}")
```

### Quality Score Visualization

```python
from dna_sequencing import fastq_parser

# Generate a quality score histogram
sequences, quality_scores = fastq_parser.readFastQ("sample.fastq")
hist, quality = fastq_parser.histogram(quality_scores)
# This saves a histogram to 'quality_histogram.png'
```

## Requirements

- Python 3.6+
- matplotlib

## License

MIT License 