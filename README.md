# SequenceForge-Lite

![Python3](https://img.shields.io/badge/Language-Python3-steelblue)
![OS](https://img.shields.io/badge/OS-_Windows_|_Mac_|_Linux-steelblue)
![License](https://img.shields.io/badge/License-MIT-steelblue)

<img src="https://github.com/iliapopov17/SequenceForge-Lite/blob/main/imgs/SeqForgeLite_logo_light.png#gh-light-mode-only" width = 50%/>
<img src="https://github.com/iliapopov17/SequenceForge-Lite/blob/main/imgs/SeqForgeLite_logo_dark.png#gh-dark-mode-only" width = 50%/>

> SequenceForge-Lite is a lightweight tool designed to work with biological sequence data, providing various functionalities for filtering FASTQ files and manipulating FASTA files. Additionally, it offers utilities for parsing BLAST output files.

## Table of contents

- [Features](#features)
- [Installation](#installation)
- [Usage Guide](#usage-guide)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Contact](#contact)

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
- Demonstrational data is available in `demo_data` folder<br>

🔗 Visit [SequenceForge-Lite tutorial](https://iliapopov17.github.io/SequenceForge-Lite/) web page

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
