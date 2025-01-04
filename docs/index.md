```markdown
# FastaFile Class Documentation

The `FastaFile` class is a Python tool for handling FASTA files. It provides methods for parsing, analyzing, and manipulating sequences in FASTA format. The class supports DNA, RNA, and protein sequences.

---

## Table of Contents

1. [Installation](#installation)
2. [Usage](#usage)
3. [API Documentation](#api-documentation)
   - [Initialization](#initialization)
   - [Sequence Analysis](#sequence-analysis)
   - [File Operations](#file-operations)
4. [Examples](#examples)
5. [Contributing](#contributing)
6. [License](#license)

---

## Installation

To use the `FastaFile` class, clone the repository and install the package:

```bash
git clone https://github.com/yourusername/fastafile-class.git
cd fastafile-class
pip install .
```

---

## Usage

Here’s a quick example of how to use the `FastaFile` class:

```python
from fastafile import FastaFile

# Initialize the FastaFile object
fasta = FastaFile('example.fasta')

# Get all entries
entries = fasta.get_entries()
for header, sequence in entries:
    print(f"Header: {header}")
    print(f"Sequence: {sequence}")

# Calculate GC content for the first sequence
gc_content = fasta.calculate_gc_content(entries[0][1])
print(f"GC Content: {gc_content:.2f}%")

# Write sequences to a new FASTA file
fasta.write_to_fasta('output.fasta')
```

---

## API Documentation

### Initialization

#### `FastaFile(file_path: str)`
- **Description**: Initializes the `FastaFile` object by reading and parsing the FASTA file.
- **Arguments**:
  - `file_path` (str): Path to the FASTA file.
- **Returns**: None.

---

### Sequence Analysis

#### `detect_sequence_type(sequence: str) -> str`
- **Description**: Detects whether the sequence is DNA, RNA, or protein.
- **Arguments**:
  - `sequence` (str): The sequence to analyze.
- **Returns**: Type of the sequence (`"DNA"`, `"RNA"`, `"Protein"`, or `"Unknown"`).

#### `calculate_gc_content(sequence: str) -> float`
- **Description**: Calculates the GC content of a DNA or RNA sequence.
- **Arguments**:
  - `sequence` (str): The sequence to analyze.
- **Returns**: GC content as a percentage.
- **Raises**: `ValueError` if the sequence is not DNA or RNA.

#### `calculate_at_content(sequence: str) -> float`
- **Description**: Calculates the AT content of a DNA or RNA sequence.
- **Arguments**:
  - `sequence` (str): The sequence to analyze.
- **Returns**: AT content as a percentage.
- **Raises**: `ValueError` if the sequence is not DNA or RNA.

#### `reverse_complement(sequence: str) -> str`
- **Description**: Generates the reverse complement of a DNA or RNA sequence.
- **Arguments**:
  - `sequence` (str): The sequence to analyze.
- **Returns**: Reverse complement of the sequence.
- **Raises**: `ValueError` if the sequence is not DNA or RNA.

#### `translate_rna_to_protein(sequence: str) -> str`
- **Description**: Translates an RNA sequence to a protein sequence.
- **Arguments**:
  - `sequence` (str): The RNA sequence to translate.
- **Returns**: Translated protein sequence.
- **Raises**: `ValueError` if the sequence is not RNA.

---

### File Operations

#### `get_entries() -> list[tuple[str, str]]`
- **Description**: Returns all entries in the FASTA file as a list of `(header, sequence)` tuples.
- **Returns**: List of tuples containing headers and sequences.

#### `write_to_fasta(output_file: str, entries: list[tuple[str, str]] = None) -> None`
- **Description**: Writes entries to a new FASTA file.
- **Arguments**:
  - `output_file` (str): Path to the output FASTA file.
  - `entries` (list[tuple[str, str]], optional): List of `(header, sequence)` tuples to write. If `None`, all entries are written.

#### `filter_sequences_by_length(min_length: int = 0, max_length: float = float('inf')) -> list[tuple[str, str]]`
- **Description**: Filters sequences based on length.
- **Arguments**:
  - `min_length` (int, optional): Minimum sequence length. Defaults to `0`.
  - `max_length` (float, optional): Maximum sequence length. Defaults to infinity.
- **Returns**: Filtered list of `(header, sequence)` tuples.

#### `find_longest_sequence() -> tuple[str, str]`
- **Description**: Finds the longest sequence in the FASTA file.
- **Returns**: Tuple containing the header and sequence of the longest entry.

#### `find_shortest_sequence() -> tuple[str, str]`
- **Description**: Finds the shortest sequence in the FASTA file.
- **Returns**: Tuple containing the header and sequence of the shortest entry.

#### `count_sequences() -> int`
- **Description**: Counts the total number of sequences in the FASTA file.
- **Returns**: Total number of sequences.

---

## Examples

### Example 1: Calculate GC Content
```python
from fastafile import FastaFile

fasta = FastaFile('example.fasta')
entries = fasta.get_entries()
gc_content = fasta.calculate_gc_content(entries[0][1])
print(f"GC Content: {gc_content:.2f}%")
```

### Example 2: Translate RNA to Protein
```python
from fastafile import FastaFile

fasta = FastaFile('example.fasta')
protein = fasta.translate_rna_to_protein("AUGGCCAUUGUAA")
print(f"Translated Protein: {protein}")
```

---

## Contributing

Contributions are welcome! Please follow these steps:
1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Submit a pull request.

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
```
