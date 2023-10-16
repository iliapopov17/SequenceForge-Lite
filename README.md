# BetaFold-3

BetaFold 3 is a versatile tool that lets you do anything you want with sequences:
1. DNA or RNA.
2. amino acids.
3. It also allows powerful filtering of fastq files.

## DNA/RNA tool

### Description

Useful tool for working with nucleic acid sequences.

The `dna_rna_tools.py` programme contains the `run_dna_rna_tools` function.
The `run_dna_rna_tools` function takes as input an arbitrary number of arguments with DNA or RNA sequences (*str*) and the name of the procedure to be executed (it is always the last argument, *str*, see example of use).
The command then performs the specified action on all submitted sequences.
If one sequence is submitted, a string with the result is returned.
If more than one sequence is submitted, a list of strings is returned. 

### List of operations:

#### 1. Transcribe DNA to RNA
The dna_rna_tool offers allows to translate a DNA sequence into an RNA sequence (transcription).
#### 2. Reversing sequence
This function allows to reverse the DNA and RNA sequence
#### 3. Complementary sequence
This function allows to display the complementary sequence of DNA or RNA
#### 4. Reverse complementary sequence
This function allows to reverse the complementary sequence of DNA or RNA

### Usage

To run dna_rna tool you need to use the function ***run_dna_rna_tools*** with the following arguments:

```python
from src.betafold_3_dna_rna_tools import run_dna_rna_tools
run_dna_rna_tools(*args: str) -> str:
```

- `*args (str):`
  1. sequence(s) - one or multiple DNA or RNA sequences
  2. procedure - DNA or RNA processing procedure

**Available procedures list**
- `transcribe` - print a transcribed sequence
- `reverse` - print reversed sequence
- `complement` - print complementary sequence
- `reverse_complement` - print reverse complementary sequence

**Usage Example**

```python
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('ATG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta'].
```

**Features:**

- The programme preserves character case (e.g. **complement AtGc** is **TaCg**)
- The programme works **only** with nucleic acid sequences. For example, the AUTGC sequence cannot exist because it contains T and U.

## Amino analyzer tool

### Description:

Useful tool for working with amino acid sequences.

### List of procedures:

#### 1. Protein molecular weight calculation
The amino_analyzer offers the capability to calculate the molecular weight of a amino acid sequence. Users can choose between average and monoisotopic weights.
#### 2. Hydrophobicity analysis
This function counts the quantity of hydrophobic and hydrophilic amino acids within a amino acid sequence.
#### 3. Cleavage site identification
Researchers can identify cleavage sites in a given peptide sequence using a specified enzyme. The tool currently supports two commonly used enzymes, trypsin and chymotrypsin.
#### 4. One-letter to three-Letter code conversion
The amino_analyzer provides a function to convert a amino acid sequence from the standard one-letter amino acid code to the three-letter code. 
#### 5. Sulphur-containing amino acid counting
The tool allows a quick determine the number of sulphur-containing amino acids, namely Cysteine (C) and Methionine (M), within a amino acid sequence. 

### Usage

To run amino_analyzer tool you need to use the function ***run_amino_analyzer*** with the following arguments:

```python
from src.betafold_3_amino_analyzer import run_amino_analyzer
run_amino_analyzer(sequence: str, procedure: str, *, weight_type: str = 'average', enzyme: str = 'trypsin'):
```

- `sequence: str:` The input amino acid sequence in one-letter code.
- `procedure: str:` The procedure to perform over your amino acid sequence.
- `weight_type: str = 'average':` default argument for `aa_weight` function. `weight_type = 'monoisotopic'` can be used as another option.
- `enzyme: str = 'trypsine':` default argument for `peptide_cutter` function. `enzyme = 'chymotrypsin'` can be used as another option
    
    
**Available procedures list**
-   `aa_weight` —  calculates the amino acids weight in a amino acid sequence.
-   `count_hydroaffinity` — counts the quantity of hydrophobic and hydrophilic amino acids in a amino acid sequence.
-   `peptide_cutter` — identifies cleavage sites in a given peptide sequence using a specified enzyme (trypsine or chymotripsine).
-   `one_to_three_letter_code` — converts a amino acid sequence from one-letter amino acid code to three-letter code.
-   `sulphur_containing_aa_counter` - counts sulphur-containing amino acids in a amino acid sequence.

You can also use each function separately by importing them in advance. Below are the available functions and their respective purposes:

##### 1. **aa_weight** function calculates the weight of amino acids in a amino acid sequence:
   The type of weight to use, either `average` or `monoisotopic`. Default is `average`.
```python
from src.betafold_3_amino_analyzer import aa_weight
aa_weight(seq: str, weight: str = `average`) -> float`
```
```python
sequence = "VLDQRKSTMA"
result = aa_weight(sequence, weight='monoisotopic')
print(result)  # Output: 1348.517
```

##### 2. **count_hydroaffinity** сounts the quantity of hydrophobic and hydrophilic amino acids in a amino acid sequence:
```python
from src.betafold_3_amino_analyzer import count_hydroaffinity 
count_hydroaffinity(seq: str) -> tuple
```
```python
sequence = "VLDQRKSTMA"
result = count_hydroaffinity(sequence)
print(result)  # Output: (3, 7)
```
##### 3. **peptide_cutter** function identifies cleavage sites in a given peptide sequence using a specified enzyme: trypsine or chymotrypsine:
```python
from src.betafold_3_amino_analyzer import peptide_cutter
peptide_cutter(sequence: str, enzyme: str = "trypsin") -> str
```
```python
sequence = "VLDQRKSTMA"
result = peptide_cutter(sequence, enzyme="trypsin")
print(result)  # Output: Found 2 trypsin cleavage sites at positions 3, 6
```
##### 4. **one_to_three_letter_code** converts a amino acid sequence from one-letter amino acid code to three-letter code.
```python
from src.betafold_3_amino_analyzer import one_to_three_letter_code
one_to_three_letter_code(sequence: str) -> str
```

```python
sequence = "VLDQRKSTMA"
result = one_to_three_letter_code(sequence)
print(result)  # Output: ValLeuAspGlnArgLysSerThrMetAla
```

##### 5. **sulphur_containing_aa_counter** counts sulphur-containing amino acids in a amino acid sequence
```python
from src.betafold_3_amino_analyzer import sulphur_containing_aa_counter
sulphur_containing_aa_counter(sequence: str) -> str
```
```python
sequence = "VLDQRKSTMA"
result = sulphur_containing_aa_counter(sequence)
print(result)  # Output: The number of sulphur-containing amino acids in the sequence is equal to 2
```

### Examples
To calculate protein molecular weight:
```python
run_amino_analyzer("VLSPADKTNVKAAW", "aa_weight")  # Output: 1481.715

run_amino_analyzer("VLSPADKTNVKAAW", "aa_weight", weight_type = 'monoisotopic')  # Output: 1480.804
```

To count hydroaffinity:
```python
run_amino_analyzer("VLSPADKTNVKAAW", "count_hydroaffinity")   # Output: (8, 6)
```

To find trypsin/chymotripsine clivage sites:
```python
run_amino_analyzer("VLSPADKTNVKAAW", "peptide_cutter") # Output: 'Found 2 trypsin cleavage sites at positions 7, 11'

run_amino_analyzer("VLSPADKTNVKAAWW", "peptide_cutter", enzyme = 'chymotrypsin') # Output: 'Found 1 chymotrypsin cleavage sites at positions 14'
```

To change to 3-letter code and count sulphur-containing amino acids.
```python
run_amino_analyzer("VLSPADKTNVKAAW", "one_to_three_letter_code") # Output: 'ValLeuSerProAlaAspLysThrAsnValLysAlaAlaTrp'

run_amino_analyzer("VLSPADKTNVKAAWM", "sulphur_containing_aa_counter") # Output: The number of sulphur-containing amino acids in the sequence is equal to 1
```

### Troubleshooting
Here are some common issues you can come ascross while using the amino-analyzer tool and their possible solutions:

1. **ValueError: Incorrect procedure**  
   If you receive this error, it means that you provided an incorrect procedure when calling `run_amino_analyzer`. Make sure you choose one of the following procedures: `aa_weight`, `count_hydroaffinity`, `peptide_cutter`, `one_to_three_letter_code`, or `sulphur_containing_aa_counter`.

   Example:
   ```python
   run_amino_analyzer("VLSPADKTNVKAAW", "incorrect_procedure")
   # Output: ValueError: Incorrect procedure. Acceptable procedures: aa_weight, count_hydroaffinity, peptide_cutter, one_to_three_letter_code, sulphur_containing_aa_counter
   ```

2. **ValueError: Incorrect sequence**
This error occurs if the input sequence provided to run_amino_analyzer contains characters that are not valid amino acids. Make sure your sequence only contains valid amino acid characters (V, I, L, E, Q, D, N, H, W, F, Y, R, K, S, T, M, A, G, P, C, v, i, l, e, q, d, n, h, w, f, y, r, k, s, t, m, a, g, p, c).

    Example:
   ```python
    run_amino_analyzer("VLSPADKTNVKAAW!", "aa_weight")
    # Output: ValueError: Incorrect sequence. Only amino acids are allowed (V, I, L, E, Q, D, N, H, W, F, Y, R, K, S, T, M, A, G, P, C, v, i, l, e, q, d, n, h,   w, f, y, r, k, s, t, m, a, g, p, c).
    ```

3. **ValueError: You have chosen an enzyme that is not provided**
This error occurs if you provide an enzyme other than "trypsin" or "chymotrypsin" when calling peptide_cutter. Make sure to use one of the specified enzymes.

    Example:
    ```python
    peptide_cutter("VLSPADKTNVKAAW", "unknown_enzyme")
    # Output: You have chosen an enzyme that is not provided. Please choose between trypsin and chymotrypsin.
    ```
4. **ValueError: You have chosen an enzyme that is not provided.**
If you encounter this error, it means that you're trying to iterate over a float value. Ensure that you're using the correct function and passing the correct arguments.

    Example:
    ```python
    result = count_hydroaffinity(123)
    # Output: TypeError: 'int' object is not iterable
    ```

## fastq tool

### Description

This utility allows you to efficiently filter DNA sequences that are written in fastq format

### List of operations:

#### 1. Check GC content
Checks if the GC content of a DNA sequence in fastq file lies within specified bounds
#### 2. Length check
Checks if the length of a DNA sequence in fastq file lies within specified bounds
#### 3. Quality check
Checks if the quality score of a DNA sequence in fastq file lies above a specified bound.

fastq_filter(seqs: dict, gc_bounds: tuple = (0, 100), length_bounds: tuple = (0, 2**32), quality_threshold: int = 0) -> dict: Filters a dictionary of sequences based on specified criteria. Returns a filtered dictionary containing only the sequences that meet the specified criteria.

### Usage

To run fastq tool you need to use the function ***run_fastq_tool*** with the following arguments:

```python
from src.betafold_3_fastq_tool import run_fastq_tool
run_fastq_tool(
    seqs: dict,
    gc_bounds: tuple = (0, 100),
    length_bounds: tuple = (0, 2**32),
    quality_threshold: int = 0,
) -> dict:
```

- `seqs: dict:` dictionary where keys are sequence identifiers and values are tuples containing sequence and quality scores
- `gc_bounds: tuple = (0, 100):` tuple with bottom and top bounds for GC content, default argument is (0,100)
- `length_bounds: tuple = (0, 2**32):` tuple with bottom and top bounds for sequence length, default argument is (0, 2**32)
- `quality_threshold: int = 0:` acceptable quality, default argument is 0

**Usage Example**

```python
EXAMPLE_FASTQ = {
    # 'name' : ('sequence', 'quality')
    '@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079804:1:SRR292678:1:1101:24563:24563': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    '@SRX079804:1:SRR292678:1:1101:30161:30161': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD'),
    '@SRX079804:1:SRR292678:1:1101:47176:47176': ('TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG', 'FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB'),
    '@SRX079804:1:SRR292678:1:1101:149302:149302': ('TAGGGTTGTATTTGCAGATCCATGGCATGCCAAAAAGAACATCGTCCCGTCCAATATCTGCAACATACCAGTTGGTTGGTA', '@;CBA=:@;@DBDCDEEE/EEEEEEF@>FBEEB=EFA>EEBD=DAEEEEB9)99>B99BC)@,@<9CDD=C,5;B::?@;A'),
    '@SRX079804:1:SRR292678:1:1101:170868:170868': ('CTGCCGAGACTGTTCTCAGACATGGAAAGCTCGATTCGCATACACTCGCTGAGTAAGAGAGTCACACCAAATCACAGATT', 'E;FFFEGFGIGGFBG;C6D<@C7CDGFEFGFHDFEHHHBBHHFDFEFBAEEEEDE@A2=DA:??C3<BCA7@DCDEG*EB'),
    '@SRX079804:1:SRR292678:1:1101:171075:171075': ('CATTATAGTAATACGGAAGATGACTTGCTGTTATCATTACAGCTCCATCGCATGAATAATTCTCTAATATAGTTGTCAT', 'HGHHHHGFHHHHFHHEHHHHFGEHFGFGGGHHEEGHHEEHBHHFGDDECEGGGEFGF<FGGIIGEBGDFFFGFFGGFGF'),
    '@SRX079804:1:SRR292678:1:1101:175500:175500': ('GACGCCGTGGCTGCACTATTTGAGGCACCTGTCCTCGAAGGGAAGTTCATCTCGACGCGTGTCACTATGACATGAATG', 'GGGGGFFCFEEEFFDGFBGGGA5DG@5DDCBDDE=GFADDFF5BE49<<<BDD?CE<A<8:59;@C.C9CECBAC=DE'),
    '@SRX079804:1:SRR292678:1:1101:190136:190136': ('GAACCTTCTTTAATTTATCTAGAGCCCAAATTTTAGTCAATCTATCAACTAAAATACCTACTGCTACTACAAGTATT', 'DACD@BEECEDE.BEDDDDD,>:@>EEBEEHEFEHHFFHH?FGBGFBBD77B;;C?FFFFGGFED.BBABBG@DBBE'),
    '@SRX079804:1:SRR292678:1:1101:190845:190845': ('CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC', 'FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE'),
    '@SRX079804:1:SRR292678:1:1101:198993:198993': ('AGTTATTTATGCATCATTCTCATGTATGAGCCAACAAGATAGTACAAGTTTTATTGCTATGAGTTCAGTACAACA', '<<<=;@B??@<>@><48876EADEG6B<A@*;398@.=BB<7:>.BB@.?+98204<:<>@?A=@EFEFFFEEFB'),
    '@SRX079804:1:SRR292678:1:1101:204480:204480': ('AGTGAGACACCCCTGAACATTCCTAGTAAGACATCTTTGAATATTACTAGTTAGCCACACTTTAAAATGACCCG', '<98;<@@@:@CD@BCCDD=DBBCEBBAAA@9???@BCDBCGF=GEGDFGDBEEEEEFFFF=EDEE=DCD@@BBC')
    }

run_fastq_tool(EXAMPLE_FASTQ,
               (20, 60),
               (0, 2**32),
               34)

# {'@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA',
  'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
 '@SRX079804:1:SRR292678:1:1101:171075:171075': ('CATTATAGTAATACGGAAGATGACTTGCTGTTATCATTACAGCTCCATCGCATGAATAATTCTCTAATATAGTTGTCAT',
  'HGHHHHGFHHHHFHHEHHHHFGEHFGFGGGHHEEGHHEEHBHHFGDDECEGGGEFGF<FGGIIGEBGDFFFGFFGGFGF'),
 '@SRX079804:1:SRR292678:1:1101:190845:190845': ('CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC',
  'FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE')}
```
