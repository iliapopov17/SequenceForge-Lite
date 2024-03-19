from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from abc import ABC, abstractmethod
from typing import Union


def fastq_filter(
    input_fastq: str = None,
    output_fastq: str = None,
    *,
    gc_bound: Union[tuple, int, float] = (0, 100),
    length_bound: Union[tuple, int, float] = (0, 2**32),
    quality_threshold: Union[int, float] = 0,
) -> None:
    """
    This function work with FASTQ files and filters them by
    GC content, length and Q-score.

    Arguments (positional):
    - input_path (str): full path to the file that you want to work with
    - output_filename (str): enter just a name of the file, don't add extention

    Arguments (keyword):
    - gc_bound (tuple, int, float): tuple of required range of GC percentage (inclusive),
    num or float if only higher border of the range is needed (exclusive).
    - length_bound (tuple, int, float): tuple of required range of sequences length (inclusive),
    num or float if only higher border of the range is needed (exclusive).
    - quality_threshold (int): int of lowest level of Q-score (inclusive).

    Output:
    - list of BioSeq records. Write file to .fastq
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
    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, index):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def check_alphabet(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    complement_dict = {"A": "", "T": "", "G": "", "C": ""}

    def __init__(self, sequence):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def check_alphabet(self):
        valid_alphabet = set("ATGC")
        return set(self.sequence) <= valid_alphabet

    def complement(self):
        complemented_sequence = [
            self.complement_dict.get(base, base) for base in self.sequence
        ]
        return self._create_instance("".join(complemented_sequence))

    def _create_instance(self, sequence):
        return self.__class__(sequence)

    def gc_content(self):
        gc_count = sum(base in "GCgc" for base in self.sequence)
        return gc_count / len(self.sequence)


class DNASequence(NucleicAcidSequence):
    complement_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def transcribe(self):
        transcription_rules = {"A": "U", "T": "A", "G": "C", "C": "G"}
        transcribed_sequence = "".join(
            transcription_rules.get(base, base) for base in self.sequence
        )
        return RNASequence(transcribed_sequence)


class RNASequence(NucleicAcidSequence):
    complement_dict = {"A": "U", "U": "A", "G": "C", "C": "G"}

    def codons(self):
        return [self.sequence[i : i + 3] for i in range(0, len(self.sequence), 3)]


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, sequence):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def check_alphabet(self):
        valid_alphabet = set("ACDEFGHIKLMNPQRSTVWY")
        return set(self.sequence) <= valid_alphabet

    def get_molecular_weight(self):
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
        return sum(molecular_weights[aa] for aa in self.sequence)
