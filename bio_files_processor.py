from typing import Optional
from dataclasses import dataclass


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: Optional[str] = None
) -> None:
    """
    Convert a multiline FASTA to a one-line FASTA.

    Parameters:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (Optional[str]): Path to the output FASTA file. Defaults to None.

    Returns:
        None
    """
    # Read the input FASTA file
    with open(input_fasta, "r") as infile:
        fasta_lines = infile.readlines()

    # Initialize variables
    sequences = {}
    current_sequence = None

    # Process each line in the FASTA file
    for line in fasta_lines:
        line = line.strip()
        if line.startswith(">"):  # Identify sequence headers
            current_sequence = line
            sequences[current_sequence] = ""
        else:  # Append lines to form a single sequence string
            sequences[current_sequence] += line

    # Set default output file name if none provided
    if output_fasta is None:
        output_fasta = input_fasta.replace(".fasta", "_oneline.fasta")

    # Write one-line FASTA sequences to the output file
    with open(output_fasta, "w") as outfile:
        for header, sequence in sequences.items():
            outfile.write(f"{header}\n")
            outfile.write(f"{sequence}\n")

    print(f"Converted multiline FASTA to one-line FASTA. Saved as {output_fasta}")


def change_fasta_start_pos(
    input_fasta: str, shift: int, output_fasta: Optional[str] = None
) -> None:
    """
    Shift the start position of all FASTA sequences in the file.

    Parameters:
        input_fasta (str): Path to the input FASTA file.
        shift (int): Positions to shift the start, can be negative.
        output_fasta (Optional[str]): Path to the output FASTA file. Defaults to None.

    Returns:
        None: The output is written to a new FASTA file.
    """
    # Create the output FASTA file name if not provided
    if output_fasta is None:
        output_fasta = input_fasta.replace(".fasta", "_shifted.fasta")

    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        header = None
        sequence = []

        for line in infile:
            line = line.strip()
            if line.startswith(">"):  # Check for header
                if header:
                    # Process the previous sequence
                    shifted_sequence = sequence[shift:] + sequence[:shift]
                    # Write the previous sequence to the output file
                    outfile.write(f"{header}\n{shifted_sequence}\n")
                # Reset for new sequence
                header = line
                sequence = ""
            else:
                # Continue reading the sequence
                sequence += line

        # Don't forget to process the last sequence in the file
        if header:
            shifted_sequence = sequence[shift:] + sequence[:shift]
            outfile.write(f"{header}\n{shifted_sequence}\n")

    print(f"Shifted FASTA sequence saved to {output_fasta}")


def parse_blast_output(input_file: str, output_file: Optional[str] = None) -> None:
    """
    Parse BLAST output to select and sort top hits.

    Parameters:
        input_file (str): Path to the input BLAST file.
        output_file (Optional[str]): Path to the output file. Defaults to None.

    Returns:
        None
    """
    # Read the input BLAST file
    with open(input_file, "r") as infile:
        lines = infile.readlines()

    # Initialize the list of best results
    best_results = []

    # Extract the best result for each query
    for i, line in enumerate(lines):
        if line.startswith("Description"):
            first_result = lines[i + 1].split("    ")[0]
            best_results.append(first_result.strip("."))

    best_results.sort()  # Sort results alphabetically

    # Set default output file name if none provided
    if output_file is None:
        output_file = input_file.replace(".txt", "_parsed.txt")

    # Write the sorted best results to the output file
    with open(output_file, "w") as outfile:
        for description in best_results:
            outfile.write(description + "\\n")


@dataclass
class FastaRecord:
    id: str
    seq: str
    description: str

    def __repr__(self) -> str:
        return f"<FastaRecord id={self.id} description='{self.description}' seq_length={len(self.seq)}>"


class OpenFasta:
    """
    Context manager for opening and iterating over a FASTA file.

    Usage:
        with OpenFasta('path_to_file.fasta') as fasta:
            for record in fasta:
                process(record)
    """

    def __init__(self, file_path: str) -> None:
        self.file_path = file_path
        self.file = None

    def __enter__(self):
        self.file = open(self.file_path, "r")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        if self.file:
            self.file.close()  # Ensure the file is closed after processing

    def __iter__(self):
        return self

    def __next__(self) -> FastaRecord:
        line = self.file.readline()
        if not line:
            raise StopIteration
        id, description = line[1:].strip().split(" ", 1)
        seq = ""
        while True:
            pos = self.file.tell()
            line = self.file.readline()
            if not line or line.startswith(">"):
                self.file.seek(pos)
                break
            seq += line.strip()
        return FastaRecord(id, seq, description)

    def read_record(self) -> FastaRecord:
        return next(self)

    def read_records(self) -> list:
        return list(self)
