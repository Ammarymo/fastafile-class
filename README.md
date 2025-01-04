# FastaFile Class

![License](https://img.shields.io/badge/license-MIT-blue.svg)

A Python class for handling FASTA files. This class provides methods for parsing, analyzing, and manipulating sequences in FASTA format. It supports DNA, RNA, and protein sequences.

---

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [License](#license)

---

## Features

- **Parse FASTA files**: Read and parse FASTA files into headers and sequences.
- **Detect sequence type**: Identify whether a sequence is DNA, RNA, or protein.
- **Calculate GC and AT content**: Compute the percentage of GC or AT bases in a sequence.
- **Reverse complement**: Generate the reverse complement of a DNA or RNA sequence.
- **Translate RNA to protein**: Convert an RNA sequence to a protein sequence using a codon table.
- **Filter sequences by length**: Extract sequences within a specified length range.
- **Find longest/shortest sequences**: Identify the longest or shortest sequence in the file.
- **Write to FASTA**: Save sequences to a new FASTA file.

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
fasta = FastaFile('examples/example.fasta')

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

For more examples, check out the [examples](examples/) directory.

---

## Documentation

For detailed documentation, including a full list of methods and their usage, see the [docs/index.md](docs/index.md) file.

---

## Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Submit a pull request.

For more information, see the [Contributing Guidelines](CONTRIBUTING.md).

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- Thanks to Me

