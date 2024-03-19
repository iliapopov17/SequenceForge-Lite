# SequenceForge-Lite

<div style='justify-content: center'>
<img src="https://github.com/iliapopov17/SequenceForge-Lite/blob/main/imgs/logo.png" align='center', width="25%">
</div>

> SequenceForge-Lite is a lightweight tool designed to work with biological sequence data, providing various functionalities for filtering FASTQ files and manipulating FASTA files. Additionally, it offers utilities for parsing BLAST output files.

## Features
### FASTQ Filtering
- Filter FASTQ files based on GC content, sequence length, and quality score.
- Specify custom ranges for GC content and sequence length.
- Set a minimum quality score threshold for sequences.
### FASTA File Manipulation
- Convert multiline FASTA files to one-line format.
- Shift the start position of one-line FASTA sequences by a specified amount.
### BLAST Output Parsing
- Extract the top hit for each query from BLAST output files.
- Results are sorted alphabetically for easy analysis.
### DNA, RNA & amino acid classes
- Calculates length of the sequence
- Calculates GC content in DNA and RNA sequences
- Prints complement sequence for DNA
- Prints transcribes RNA sequence for DNA
- Prints RNA sequence in codons
- Calculates molecular weight of amino acid sequence

## Installation

```bash
git clone https://github.com/iliapopov17/SequenceForge-Lite.git && cd SequenceForge-Lite
```

```bash
pip install biopython
```

## Usage Guide
- Demonstrational python notebook is available in `demo.ipynb` file
- Demonstrational data is available in `demo_data` folder

### Example of using the `convert_multiline_fasta_to_online` function

```python
input_fasta_file = "demo_data/example_multiline_fasta.fasta"
output_fasta_file = "demo_data/example_oneline_fasta.fasta"
convert_multiline_fasta_to_oneline(input_fasta_file, output_fasta_file)
```

```
Converted multiline FASTA to one-line FASTA. Saved as demo_data/example_oneline_fasta.fasta
```

### Example of using the `change_fasta_start_pos` function

```python
input_fasta_file = "demo_data/example_oneline_fasta.fasta"
shift_amount = 10
change_fasta_start_pos(input_fasta_file, shift_amount)
```

```
Shifted FASTA sequence saved to demo_data/example_oneline_fasta_shifted.fasta
```

### Example of using the `parse_blast_output` function

```python
input_file = "demo_data/example_blast_results.txt"
parse_blast_output(input_file)
```

### **Example of using the `fastq_filter` function**

```python
input_file = "demo_data/example_fastq.fastq"
fastq_filter(input_file, gc_bound=(40,60), length_bound=(0, 200), quality_threshold=25)
```

```
Filtered FastQ. Saved as demo_data/example_fastq_filtered.fastq
```

### Example usage of `DNASequence` class

```python
dna_sequence = DNASequence("ATCGATCG")
print("DNA Sequence:", dna_sequence)
print("Length:", len(dna_sequence))
print("GC Content:", dna_sequence.gc_content())
print("Complement:", dna_sequence.complement())
print("Transcribed RNA Sequence:", dna_sequence.transcribe())
```

```
DNA Sequence: ATCGATCG
Length: 8
GC Content: 0.5
Complement: TAGCTAGC
Transcribed RNA Sequence: UAGCUAGC
```

### Example usage of `RNASequence` class

```python
rna_sequence = RNASequence("AUCGAUCGA")
print("RNA Sequence:", rna_sequence)
print("Length:", len(rna_sequence))
print("GC Content:", rna_sequence.gc_content())
print("Codons:", rna_sequence.codons())
```

```
RNA Sequence: AUCGAUCGA
Length: 9
GC Content: 0.4444444444444444
Codons: ['AUC', 'GAU', 'CGA']
```

### Example usage of `AminoAcidSequence` class

```python
amino_acid_sequence = AminoAcidSequence("ACDEFG")
print("Amino Acid Sequence:", amino_acid_sequence)
print("Length:", len(amino_acid_sequence))
print("Molecular Weight:", amino_acid_sequence.get_molecular_weight())
```

```
Amino Acid Sequence: ACDEFG
Length: 6
Molecular Weight: 730.74
```

## Troubleshooting:
Common Issues and Solutions:
1. **File Not Found Error:**
- **Issue**: The script raises a FileNotFoundError when trying to access the input file.
- **Solution**: Verify that the input file path provided to the function is correct and the file exists in the specified location.
2. **Incorrect File Format**:
- **Issue**: The function fails to process the file due to incorrect formatting.
- **Solution**: Ensure that the input files are properly formatted according to the specifications mentioned in the function or class descriptions.

## Contributing
Contributions are welcome! If you have any ideas, bug fixes, or enhancements, feel free to open an issue or submit a pull request.

## Contact
For any inquiries or support, feel free to contact me via [email](iljapopov@gmail.com).

Happy sequencing! 🧬🔬
