# SequenceForge-Lite

![Python3](https://img.shields.io/badge/Language-Python3-steelblue)
![OS](https://img.shields.io/badge/OS-_Windows_|_Mac_|_Linux-steelblue)
![License](https://img.shields.io/badge/License-MIT-steelblue)

<div style='justify-content: center'>
<img src="https://github.com/iliapopov17/SequenceForge-Lite/blob/main/imgs/logo.png" align='center', width="25%">
</div>

> SequenceForge-Lite is a lightweight tool designed to work with biological sequence data, providing various functionalities for filtering FASTQ files and manipulating FASTA files. Additionally, it offers utilities for parsing BLAST output files.

## Table of contents

- [Features](#features)
- [Installation](#installation)
- [Usage Guide](#usage-guide)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)

## Features
### FASTQ Filtering
- Filter FASTQ files based on GC content, sequence length, and quality score.
- Specify custom ranges for GC content and sequence length.
- Set a minimum quality score threshold for sequences.
### FASTA File Manipulation
- Get quick info on each sequence in FASTA file.
- Convert multiline FASTA files to one-line format.
- Shift the start position of one-line FASTA sequences by a specified amount.
### BLAST Output Parsing
- Extract the top hit for each query from BLAST output files.
- Results are sorted alphabetically for easy analysis.
### DNA, RNA & amino acid classes
- Calculates GC content in DNA and RNA sequences
- Prints complement sequence for DNA
- Transcribes DNA sequence to RNA
- Prints RNA sequence in codons
- Finds motifs in nucleic acids sequences
- Translates RNA sequence to amino acid (without biological meaning, it does it "dumbly")
- Calculates molecular weight of amino acid sequence
### Custom RandomForestClassifier
- Self written implementation of RandomForestClassifier
- Has parallelisation functionality (speeds up 2 times when specifying 2 threads)

## Installation

```bash
git clone https://github.com/iliapopov17/SequenceForge-Lite.git && cd SequenceForge-Lite
```

```bash
pip install -r requirements.txt
```

## Usage Guide
- Demonstrational python notebook is available in `demo.ipynb` file
- Demonstrational data is available in `demo_data` folder

> DISCLAMER: I believe 'Usage Guide' must be in the README! I have read HW description, I know that it is not welcome to perserve it in README. Yet my opinion is my opinion and it cannot be changed. I plan to master my skills in writing detailed instructions and oneday I will deploy github pages with the instructions how to use this tool. But for now 'Usage Guide' will stay in README. User must be treated with love.

### Import modules

**_Input_**

```python
from bio_files_processor import *
from sequence_forge import *
```

### Example of using the `convert_multiline_fasta_to_online` function

**_Input_**

```python
! head -6 demo_data/example_multiline_fasta.fasta
```

**_Output_**

```
>GTD323452 5S_rRNA NODE_272_length_223_cov_0.720238:18-129(+)
ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCT
GGTGCTGTG
>GTD678345 16S_rRNA NODE_80_length_720_cov_1.094737:313-719(+)
TTGGCTTCTTAGAGGGACTTTTGATGTTTAATCAAAGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACT
GAGCCCTTGGGAGTGGTCCATTTGAGCCGGCAACGGCACGTTTGGACTGCAAACTTGGGCAAACTTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGT
```

**_Input_**

```python
input_fasta_file = "demo_data/example_multiline_fasta.fasta"
output_fasta_file = "demo_data/example_oneline_fasta.fasta"
convert_multiline_fasta_to_oneline(input_fasta_file, output_fasta_file)
```

**_Output_**

```
Converted multiline FASTA to one-line FASTA. Saved as demo_data/example_oneline_fasta.fasta
```

**_Input_**

```python
! head -4 demo_data/example_oneline_fasta.fasta
```

**_Output_**

```
>GTD323452 5S_rRNA NODE_272_length_223_cov_0.720238:18-129(+)
ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCTGGTGCTGTG
>GTD678345 16S_rRNA NODE_80_length_720_cov_1.094737:313-719(+)
TTGGCTTCTTAGAGGGACTTTTGATGTTTAATCAAAGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGAGCCCTTGGGAGTGGTCCATTTGAGCCGGCAACGGCACGTTTGGACTGCAAACTTGGGCAAACTTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGT
```

### Example of using the `change_fasta_start_pos` function

**_Input_**

```python
input_fasta_file = "demo_data/example_oneline_fasta.fasta"
shift_amount = 10
change_fasta_start_pos(input_fasta_file, shift_amount)
```

**_Output_**

```
Shifted FASTA sequence saved to demo_data/example_oneline_fasta_shifted.fasta
```

**_Input_**

```python
! head -4 demo_data/example_oneline_fasta_shifted.fasta
```

**_Output_**

```
>GTD323452 5S_rRNA NODE_272_length_223_cov_0.720238:18-129(+)
GACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCTGGTGCTGTGACGGCCATAG
>GTD678345 16S_rRNA NODE_80_length_720_cov_1.094737:313-719(+)
AGAGGGACTTTTGATGTTTAATCAAAGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGAGCCCTTGGGAGTGGTCCATTTGAGCCGGCAACGGCACGTTTGGACTGCAAACTTGGGCAAACTTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTGGCTTCTT
```

### Example of using the `parse_blast_output` function

**_Input_**

```python
input_file = "demo_data/example_blast_results.txt"
parse_blast_output(input_file)
```

**_Output_**

```
Best BLAST results saved to demo_data/example_blast_results_parsed.txt
```

**_Input_**

```python
! head -4 demo_data/example_blast_results_parsed.txt
```

**_Output_**

```
DNA methylase [Enterobacteriaceae]
DUF1380 domain-containing protein [Escherichia coli]
DUF1380 family protein [Enterobacteriaceae]
DUF4158 domain-containing protein [Klebsiella pneumoniae]
```

### Example of using the `fastq_filter` function

**_Input_**

```python
! head -4 demo_data/example_fastq.fastq
```

**_Output_**

```
@SRX079804:1:SRR292678:1:1101:24563:24563 1:N:0:1 BH:failed
ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG
+
BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D
```

**_Input_**

```python
input_file = "demo_data/example_fastq.fastq"
fastq_filter(input_file, gc_bound=(40,60), length_bound=(0, 50), quality_threshold=30)
```

**_Output_**

```
Filtered FastQ. Saved as demo_data/example_fastq_filtered.fastq
```

**_Input_**

```python
! head -4 demo_data/example_fastq_filtered.fastq
```

**_Output_**

```
@SRX079804:1:SRR292678:1:1101:654270:654270 1:N:0:1 BH:failed
ATCTTTCTCTCTCTGCTGCATTCTCCGCTTCAGCTCCTCAATTTCAATCA
+
GGGGGGEDGGFEGDGCGGBCEFGGGGFBFGGEGGGGGEGFCFEEGGFEDB
```

### Example usage of `DNASequence` class

**_Input_**

```python
dna_sequence = DNASequence("ACCGGCTAATCGGCT")
motif_to_find = "CGG"
print(type(dna_sequence))
print("DNA Sequence:", dna_sequence)
print("Length:", len(dna_sequence))
print("GC Content:", dna_sequence.gc_content())
print("Complement:", dna_sequence.complement())
print("Transcribed RNA Sequence:", dna_sequence.transcribe())
print(f"Indexes of {motif_to_find} motif occurrences:", dna_sequence.find_motif(motif_to_find))
```

**_Output_**

```
<class 'sequence_forge.DNASequence'>
DNA Sequence: ACCGGCTAATCGGCT
Length: 15
GC Content: 0.6
Complement: TGGCCGATTAGCCGA
Transcribed RNA Sequence: UGGCCGAUUAGCCGA
Indexes of CGG motif occurrences: [2, 10]
```

### Example usage of `RNASequence` class

**_Input_**

```python
rna_sequence = dna_sequence.transcribe()
motif_to_find = "GCC"
print(type(rna_sequence))
print("RNA Sequence:", rna_sequence)
print("Length:", len(rna_sequence))
print("GC Content:", rna_sequence.gc_content())
print("Codons:", rna_sequence.codons())
print(f"Indexes of {motif_to_find} motif occurrences:", rna_sequence.find_motif(motif_to_find))
print("Tranlated to Amino Acid Sequence:", rna_sequence.translate())
```

**_Output_**

```
<class 'sequence_forge.RNASequence'>
RNA Sequence: UGGCCGAUUAGCCGA
Length: 15
GC Content: 0.6
Codons: ['UGG', 'CCG', 'AUU', 'AGC', 'CGA']
Indexes of GCC motif occurrences: [2, 10]
Tranlated to Amino Acid Sequence: WPISR
```

### Example usage of `AminoAcidSequence` class

**_Input_**

```python
amino_acid_sequence = rna_sequence.translate()
print(type(amino_acid_sequence))
print("Amino Acid Sequence:", amino_acid_sequence)
print("Length:", len(amino_acid_sequence))
print("Molecular Weight:", amino_acid_sequence.get_molecular_weight())
```

**_Output_**

```
<class 'sequence_forge.AminoAcidSequence'>
Amino Acid Sequence: WPISR
Length: 5
Molecular Weight: 729.8299999999999
```

### Testing RandomForestClassifierCustom

Import needed library and module

**_Input_**

```python
from custom_random_forest import RandomForestClassifierCustom
from sklearn.datasets import make_classification
```

Create dataset

**_Input_**

```python
X, y = make_classification(n_samples=100000)
rf = RandomForestClassifierCustom(n_estimators=10, max_depth=10, max_features=5, random_state=42)
```

Fit using 2 `n_jobs`

**_Input_**

```python
%%time
rf.fit(X, y, n_jobs=2)
```

**_Output_**

```
CPU times: total: 109 ms
Wall time: 7.95 s
```

Predict using 2 `n_jobs`

**_Input_**

```python
%%time
predictions_2 = rf.predict_proba(X, n_jobs=2)
```

**_Output_**

```
CPU times: total: 78.1 ms
Wall time: 2.3 s
```

Fit using 1 `n_jobs`

**_Input_**

```python
%%time
rf.fit(X, y, n_jobs=1)
```

**_Output_**

```
CPU times: total: 9.56 s
Wall time: 9.57 s
```

Predict using 1 `n_jobs`

**_Input_**

```python
%%time
predictions_2 = rf.predict_proba(X, n_jobs=1)
```

**_Output_**

```
CPU times: total: 172 ms
Wall time: 171 ms
```

Time comparison (CPU times)

Time comparison (CPU times)
| |fit|predict|
|-|-|-|
|n_jobs=1|9.44 s|172 ms|
|n_jobs=2|109 ms|78.1 ms|

Finally, check that the resulting predictions match

**_Input_**

```python
predictions_1 == predictions_2
```

**_Output_**

```
array([[ True,  True],
       [ True,  True],
       [ True,  True],
       ...,
       [ True,  True],
       [ True,  True],
       [ True,  True]])
```

## Troubleshooting
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
For any inquiries or support, feel free to contact me via [email](mailto:iljapopov17@gmail.com)

Happy sequencing! 🧬🔬
