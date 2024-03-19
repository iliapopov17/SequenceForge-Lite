def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None):
    """
    Converts a multiline FASTA file to a one-line FASTA file.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str, optional): Path to the output FASTA file. If not provided, the output file will have the same name as the input file with the extension .fasta.

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
        if line.startswith(">"):
            # Save the header
            current_sequence = line
            sequences[current_sequence] = ""
        else:
            # Concatenate the sequence
            sequences[current_sequence] += line

    # Create the output FASTA file name
    if output_fasta is None:
        output_fasta = input_fasta.replace(".fasta", "_oneline.fasta")

    # Write one-line FASTA sequences to the output file
    with open(output_fasta, "w") as outfile:
        for header, sequence in sequences.items():
            outfile.write(f"{header}\n")
            outfile.write(f"{sequence}\n")

    print(f"Converted multiline FASTA to one-line FASTA. Saved as {output_fasta}")


def change_fasta_start_pos(input_fasta, shift, output_fasta=None):
    """
    Shifts the start position of a one-line FASTA sequence by the specified amount.

    Args:
        input_fasta (str): Path to the input FASTA file.
        shift (int): Number of positions to shift the start position (positive or negative).
        output_fasta (str, optional): Path to the output FASTA file. If not provided, the output file will have the same name as the input file with the extension .fasta.

    Returns:
        None
    """
    # Read the input FASTA file
    with open(input_fasta, "r") as infile:
        fasta_lines = infile.readlines()

    # Extract the sequence
    sequence = fasta_lines[1].strip()

    # Shift the sequence
    shifted_sequence = sequence[shift:] + sequence[:shift]

    # Create the output FASTA file name
    if output_fasta is None:
        output_fasta = input_fasta.replace(".fasta", "_shifted.fasta")

    # Write the shifted sequence to the output FASTA file
    with open(output_fasta, "w") as outfile:
        outfile.write(f">{fasta_lines[0].strip()}\n")
        outfile.write(f"{shifted_sequence}\n")

    print(f"Shifted FASTA sequence saved to {output_fasta}")


def parse_blast_output(input_file: str, output_file: str = None):
    """
    Selects the top hit for each query from the input file, saves them to the output file sorted alphabetically.

    Args:
        input_file (str): Path to the input file.
        output_file (str, optional): Path to the output file. If not provided, the results will be printed to the console.

    Returns:
        None
    """
    # Read the input file
    with open(input_file, "r") as infile:
        lines = infile.readlines()

    # Initialize the list of best results
    best_results = []

    # Extract the best result for each query
    for i, line in enumerate(lines):
        if line.startswith("Description"):
            first_result = lines[i + 1].split("    ")[0]
            best_results.append(first_result.strip("."))

    # Sort the results
    best_results.sort()

    # Create the output file name
    if output_file is None:
        output_file = input_file.replace(".txt", "_parsed.txt")

    # Write the best results to the output file or print to console
    with open(output_file, "w") as outfile:
        for description in best_results:
            outfile.write(description + "\n")

    print(f"Best BLAST results saved to {output_file}")
