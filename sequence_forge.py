from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from abc import ABC, abstractmethod
from typing import Union, List, Type


def fastq_filter(
    input_fastq: str = None,
    output_fastq: str = None,
    *,
    gc_bound: Union[tuple, int, float] = (0, 100),
    length_bound: Union[tuple, int, float] = (0, 2**32),
    quality_threshold: Union[int, float] = 0,
) -> None:
    """
    Filter FASTQ records based on GC content, sequence length, and quality score.

    Parameters:
        input_fastq (str): Path to the input FASTQ file.
        output_fastq (str): Path to the output FASTQ file.
        gc_bound (Union[Tuple[int, int], int]): GC content bounds (inclusive) or minimum GC percentage.
        length_bound (Union[Tuple[int, int], int]): Sequence length bounds (inclusive) or minimum length.
        quality_threshold (int): Minimum quality score (inclusive).

    Returns:
        None: Writes filtered sequences to an output FASTQ file.
    """

    # Create the output FASTA file name
    if output_fastq is None:
        output_fastq = input_fastq.replace(".fastq", "_filtered.fastq")

    # Create dict from FASTQ
    seqs = list(SeqIO.parse(input_fastq, "fastq"))

    # Create filtered list
    filtered_fastq = []

    for line in seqs:
        if isinstance(gc_bound, (int, float)):
            gc_check = gc_fraction(line) * 100 >= gc_bound
        else:
            gc_check = (
                gc_fraction(line) * 100 < gc_bound[0]
                or gc_fraction(line) * 100 > gc_bound[1]
            )
        if isinstance(length_bound, (int, float)):
            len_check = len(line) >= length_bound
        else:
            len_check = len(line) < length_bound[0] or len(line) > length_bound[1]
        quality_check = (
            sum(line.letter_annotations["phred_quality"])
            / len(line.letter_annotations["phred_quality"])
            < quality_threshold
        )
        if not (gc_check or len_check or quality_check):
            filtered_fastq.append(line)

    # Write  filtered data into new .fastq file
    SeqIO.write(filtered_fastq, output_fastq, "fastq")

    print(f"Filtered FastQ. Saved as {output_fastq}")


class BiologicalSequence(ABC):
    """
    Abstract base class for biological sequences, enforcing methods for subclasses.
    """

    @abstractmethod
    def __len__(self) -> int:
        """
        Return the length of the sequence.

        Returns:
            int: The number of characters in the sequence.
        """
        pass

    @abstractmethod
    def __getitem__(self, index: int):
        """
        Retrieve an element from the sequence by index.

        Parameters:
            index (int): The position in the sequence to retrieve.

        Returns:
            The element at the specified index in the sequence.
        """
        pass

    @abstractmethod
    def __str__(self) -> str:
        """
        Return a string representation of the sequence.

        Returns:
            str: The sequence as a standard string.
        """
        pass

    @abstractmethod
    def check_alphabet(self) -> bool:
        """
        Check if the sequence contains only valid characters.

        Returns:
            bool: True if all characters in the sequence are valid, False otherwise.
        """
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Base class for nucleic acid sequences, supporting basic operations like obtaining a complement
    and calculating GC content.
    """

    complement_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
    }

    def __init__(self, sequence: str) -> None:
        """
        Initialize a NucleicAcidSequence with a sequence string.

        Parameters:
            sequence (str): The nucleotide sequence.
        """
        self.sequence = sequence

    def __len__(self) -> int:
        """
        Return the length of the sequence.

        Returns:
            int: The number of characters in the sequence.
        """
        return len(self.sequence)

    def __getitem__(self, index: int) -> str:
        """
        Retrieve a nucleotide from the sequence by index.

        Parameters:
            index (int): The position in the sequence to retrieve.

        Returns:
            str: The nucleotide at the specified index.
        """
        return self.sequence[index]

    def __str__(self) -> str:
        """
        Return a string representation of the nucleic acid sequence.

        Returns:
            str: The nucleotide sequence as a standard string.
        """
        return self.sequence

    def check_alphabet(self) -> bool:
        """
        Check if the sequence contains only valid nucleotides (A, T, G, C).

        Returns:
            bool: True if all characters in the sequence are valid, False otherwise.
        """
        valid_alphabet = set("AUTGC")
        return set(self.sequence) <= valid_alphabet

    def complement(self) -> "NucleicAcidSequence":
        """
        Compute the complement of the nucleic acid sequence.

        Returns:
            NucleicAcidSequence: A new instance of NucleicAcidSequence representing the complement.
        """
        complemented_sequence = [
            self.complement_dict.get(base, base) for base in self.sequence
        ]
        return self._create_instance("".join(complemented_sequence))

    def _create_instance(self, sequence: str) -> "NucleicAcidSequence":
        """
        Create a new instance of the same class with a given sequence.

        Parameters:
            sequence (str): The nucleotide sequence for the new instance.

        Returns:
            NucleicAcidSequence: A new instance of this class.
        """
        return self.__class__(sequence)

    def gc_content(self) -> float:
        """
        Calculate the GC content of the sequence.

        Returns:
            float: The percentage of nucleotides in the sequence that are G or C.
        """
        gc_count = sum(1 for base in self.sequence if base in "GC")
        return (gc_count / len(self.sequence)) * 100 if self.sequence else 0

    def find_motif(self, motif: str) -> List[int]:
        """
        Find all occurrences of a motif in the sequence.

        Parameters:
            motif (str): The motif to search for.

        Returns:
            List[int]: A list of start indices where the motif appears in the sequence.
        """
        motif_indices = []
        motif_length = len(motif)
        sequence_length = len(self)

        for i in range(sequence_length - motif_length + 1):
            if self.sequence[i : i + motif_length] == motif:
                motif_indices.append(i)

        return motif_indices


class DNASequence(NucleicAcidSequence):
    """
    Represents a DNA sequence and provides DNA-specific operations like transcription.
    """

    complement_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def transcribe(self) -> "RNASequence":
        """
        Transcribe the DNA sequence to RNA.

        This method converts the DNA sequence into its RNA equivalent based on standard
        DNA-to-RNA transcription rules.

        Returns:
            RNASequence: The transcribed RNA sequence.
        """
        transcription_rules = {"A": "U", "T": "A", "G": "C", "C": "G"}
        transcribed_sequence = "".join(
            transcription_rules.get(base, base) for base in self.sequence
        )
        return RNASequence(transcribed_sequence)


class RNASequence(NucleicAcidSequence):
    """
    Represents an RNA sequence and provides methods for processing RNA into proteins.
    """

    complement_dict = {"A": "U", "U": "A", "G": "C", "C": "G"}

    codon_table = {
        "UUU": "F",
        "UUC": "F",
        "UUA": "L",
        "UUG": "L",
        "CUU": "L",
        "CUC": "L",
        "CUA": "L",
        "CUG": "L",
        "AUU": "I",
        "AUC": "I",
        "AUA": "I",
        "AUG": "M",
        "GUU": "V",
        "GUC": "V",
        "GUA": "V",
        "GUG": "V",
        "UCU": "S",
        "UCC": "S",
        "UCA": "S",
        "UCG": "S",
        "CCU": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "ACU": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "GCU": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "UAU": "Y",
        "UAC": "Y",
        "UAA": "*",
        "UAG": "*",
        "CAU": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "AAU": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "GAU": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "UGU": "C",
        "UGC": "C",
        "UGA": "*",
        "UGG": "W",
        "CGU": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "AGU": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GGU": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }

    def codons(self) -> list:
        """
        Generates a list of codons from the RNA sequence.

        Returns:
            list: List of codon strings.
        """
        return [self.sequence[i : i + 3] for i in range(0, len(self.sequence), 3)]

    def translate(self) -> "AminoAcidSequence":
        """
        Translate the RNA sequence into an amino acid sequence using the genetic code.

        Returns:
            AminoAcidSequence: Represents the translated amino acid sequence.
        """
        codons = [self.sequence[i : i + 3] for i in range(0, len(self.sequence), 3)]
        translated_sequence = ""
        for codon in codons:
            if len(codon) == 3:
                translated_sequence += RNASequence.codon_table.get(codon, "?")
        return AminoAcidSequence(translated_sequence)


class AminoAcidSequence(BiologicalSequence):
    """
    Represents an amino acid sequence and provides methods for biochemical properties.
    """

    def __init__(self, sequence: str) -> None:
        """
        Initialize the AminoAcidSequence with a sequence string.

        Parameters:
            sequence (str): The amino acid sequence.
        """
        self.sequence = sequence

    def __len__(self) -> int:
        """
        Get the length of the amino acid sequence.

        Returns:
            int: Length of the sequence.
        """
        return len(self.sequence)

    def __getitem__(self, index: int) -> str:
        """
        Get the amino acid at a specific index.

        Parameters:
            index (int): Index of the amino acid in the sequence.

        Returns:
            str: Amino acid at the specified index.
        """
        return self.sequence[index]

    def __str__(self) -> str:
        """
        Get the string representation of the amino acid sequence.

        Returns:
            str: The amino acid sequence.
        """
        return self.sequence

    def check_alphabet(self) -> bool:
        """
        Check if the sequence contains only valid amino acids.

        Returns:
            bool: True if the sequence contains only valid amino acids, False otherwise.
        """
        valid_alphabet = set("ACDEFGHIKLMNPQRSTVWY")
        return set(self.sequence) <= valid_alphabet

    def get_molecular_weight(self) -> float:
        """
        Calculate the total molecular weight of the amino acid sequence.

        Returns:
            float: Total molecular weight of the sequence.
        """
        # Amino acid to molecular weight mapping (abbreviated for space)
        molecular_weights = {
            "A": 89.09,
            "C": 121.16,
            "D": 133.10,
            "E": 147.13,
            "F": 165.19,
            "G": 75.07,
            "H": 155.16,
            "I": 131.18,
            "K": 146.19,
            "L": 131.18,
            "M": 149.21,
            "N": 132.12,
            "P": 115.13,
            "Q": 146.15,
            "R": 174.20,
            "S": 105.09,
            "T": 119.12,
            "V": 117.15,
            "W": 204.23,
            "Y": 181.19,
        }
        # Calculate weight by summing the weights of each amino acid in the sequence
        return sum(molecular_weights[aa] for aa in self.sequence)
