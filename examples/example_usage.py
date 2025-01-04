from fastafile.fastafile import FastaFile

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
