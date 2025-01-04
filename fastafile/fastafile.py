class FastaFile:
    def __init__(self, file_path: str) -> None:
        """
        Initialize the FastaFile object by reading and parsing the FASTA file.

        Args:
            file_path (str): Path to the FASTA file.
        """
        self.file_path = file_path
        self.entries = self._parse_fasta()

    def _read_fasta(self) -> str:
        """
        Read the FASTA file and return its content as a string.

        Returns:
            str: Content of the FASTA file.
        """
        with open(self.file_path, 'r') as file:
            return file.read()

    def _parse_fasta(self) -> list[tuple[str, str]]:
        """
        Parse the FASTA file content into a list of (header, sequence) tuples.

        Returns:
            list[tuple[str, str]]: List of tuples containing headers and sequences.
        """
        file_content = self._read_fasta()
        entries = file_content.split('>')[1:]  # Skip the first element (it's empty)
        fasta_data = []
        for entry in entries:
            lines = entry.split('\n')
            header = lines[0].strip()
            sequence = ''.join(lines[1:]).strip()
            fasta_data.append((header, sequence))
        return fasta_data

    def detect_sequence_type(self, sequence: str) -> str:
        """
        Detect whether the sequence is DNA, RNA, or protein.

        Args:
            sequence (str): The sequence to analyze.

        Returns:
            str: Type of the sequence ("DNA", "RNA", "Protein", or "Unknown").
        """
        sequence = sequence.upper()
        dna_bases = {'A', 'T', 'C', 'G'}
        rna_bases = {'A', 'U', 'C', 'G'}
        protein_bases = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}

        if all(base in dna_bases for base in sequence):
            return "DNA"
        elif all(base in rna_bases for base in sequence):
            return "RNA"
        elif all(base in protein_bases for base in sequence):
            return "Protein"
        else:
            return "Unknown"

    def calculate_gc_content(self, sequence: str) -> float:
        """
        Calculate the GC content of a DNA or RNA sequence.

        Args:
            sequence (str): The sequence to analyze.

        Returns:
            float: GC content as a percentage.

        Raises:
            ValueError: If the sequence is not DNA or RNA.
        """
        sequence = sequence.upper()
        sequence_type = self.detect_sequence_type(sequence)
        
        if sequence_type == "DNA":
            gc_count = sequence.count('G') + sequence.count('C')
        elif sequence_type == "RNA":
            gc_count = sequence.count('G') + sequence.count('C')
        else:
            raise ValueError("GC content is only applicable to DNA or RNA sequences.")
        
        total_length = len(sequence)
        if total_length == 0:
            return 0.0
        gc_content = (gc_count / total_length) * 100
        return gc_content

    def calculate_at_content(self, sequence: str) -> float:
        """
        Calculate the AT content of a DNA or RNA sequence.

        Args:
            sequence (str): The sequence to analyze.

        Returns:
            float: AT content as a percentage.

        Raises:
            ValueError: If the sequence is not DNA or RNA.
        """
        sequence = sequence.upper()
        sequence_type = self.detect_sequence_type(sequence)
        
        if sequence_type == "DNA":
            at_count = sequence.count('A') + sequence.count('T')
        elif sequence_type == "RNA":
            at_count = sequence.count('A') + sequence.count('U')
        else:
            raise ValueError("AT content is only applicable to DNA or RNA sequences.")
        
        total_length = len(sequence)
        if total_length == 0:
            return 0.0
        at_content = (at_count / total_length) * 100
        return at_content

    def reverse_complement(self, sequence: str) -> str:
        """
        Generate the reverse complement of a DNA or RNA sequence.

        Args:
            sequence (str): The sequence to analyze.

        Returns:
            str: Reverse complement of the sequence.

        Raises:
            ValueError: If the sequence is not DNA or RNA.
        """
        sequence = sequence.upper()
        sequence_type = self.detect_sequence_type(sequence)
        
        if sequence_type == "DNA":
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        elif sequence_type == "RNA":
            complement = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
        else:
            raise ValueError("Reverse complement is only applicable to DNA or RNA sequences.")
        
        return ''.join([complement[base] for base in reversed(sequence)])

    def translate_rna_to_protein(self, sequence: str) -> str:
        """
        Translate an RNA sequence to a protein sequence.

        Args:
            sequence (str): The RNA sequence to translate.

        Returns:
            str: Translated protein sequence.

        Raises:
            ValueError: If the sequence is not RNA.
        """
        sequence = sequence.upper()
        if self.detect_sequence_type(sequence) != "RNA":
            raise ValueError("Sequence must be RNA to translate to protein.")
        
        codon_table = {
            'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
            'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
            'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
            'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
            'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
            'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
            'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
            'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
            'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
            'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W'
        }
        protein = []
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            protein.append(codon_table.get(codon, '?'))  # '?' for unknown codons
        return ''.join(protein)

    def get_entries(self) -> list[tuple[str, str]]:
        """
        Return all entries in the FASTA file as a list of (header, sequence) tuples.

        Returns:
            list[tuple[str, str]]: List of tuples containing headers and sequences.
        """
        return self.entries

    def write_to_fasta(self, output_file: str, entries: list[tuple[str, str]] = None) -> None:
        """
        Write entries to a new FASTA file.

        Args:
            output_file (str): Path to the output FASTA file.
            entries (list[tuple[str, str]], optional): List of (header, sequence) tuples to write.
                If None, all entries are written. Defaults to None.
        """
        if entries is None:
            entries = self.entries
        with open(output_file, 'w') as file:
            for header, sequence in entries:
                file.write(f">{header}\n{sequence}\n")

    def filter_sequences_by_length(self, min_length: int = 0, max_length: float = float('inf')) -> list[tuple[str, str]]:
        """
        Filter sequences based on length.

        Args:
            min_length (int, optional): Minimum sequence length. Defaults to 0.
            max_length (float, optional): Maximum sequence length. Defaults to infinity.

        Returns:
            list[tuple[str, str]]: Filtered list of (header, sequence) tuples.
        """
        filtered_entries = []
        for header, sequence in self.entries:
            if min_length <= len(sequence) <= max_length:
                filtered_entries.append((header, sequence))
        return filtered_entries

    def find_longest_sequence(self) -> tuple[str, str]:
        """
        Find the longest sequence in the FASTA file.

        Returns:
            tuple[str, str]: Tuple containing the header and sequence of the longest entry.
        """
        return max(self.entries, key=lambda x: len(x[1]))

    def find_shortest_sequence(self) -> tuple[str, str]:
        """
        Find the shortest sequence in the FASTA file.

        Returns:
            tuple[str, str]: Tuple containing the header and sequence of the shortest entry.
        """
        return min(self.entries, key=lambda x: len(x[1]))

    def count_sequences(self) -> int:
        """
        Count the total number of sequences in the FASTA file.

        Returns:
            int: Total number of sequences.
        """
        return len(self.entries)
    
