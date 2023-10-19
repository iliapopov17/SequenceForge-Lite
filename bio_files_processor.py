import os

def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> None:
    '''
    This function converts multiline fasta to oneline fasta

    Parameters:
    - input_fasta (str): path to the intact multiline fasta
    - output_fasta (str): name of the output oneline fasta file
    '''

    if output_fasta is None:
        output_fasta = os.path.basename(input_fasta) + '.fasta'
    with open(input_fasta) as file_input, open(output_fasta, 'w') as file_output:
        seq = []
        for line in file_input:
            if line.startswith('>'):
                if seq:
                    file_output.write(''.join(seq) + '\n')
                seq = []
                file_output.write(line)
            else:
                seq.append(line.strip())
        file_output.write(''.join(seq))

def get_gene_sequence(file_input):
    '''
    This function parses a GenBank (GBK) file and extracts information about a gene, specifically the gene's name and its associated amino acid sequence.
    It does so by reading through the lines of the input GenBank file and identifying the "/gene" and "/translation" sections, capturing the gene name and sequence accordingly

    Parameters:
    - file_input: file-like object that provides access to the GenBank file. This object should be opened and ready for reading

    Returns:
    - gene: name of the gene, which is extracted from the GenBank file
    - seq: amino acid sequence associated with the gene, which is also extracted from the GenBank file
    '''

    gene = None
    seq = None
    for line in file_input:
        line = line.strip()
        if line.startswith('/gene'):
            gene = line[7:-1]
        elif line.startswith('/translation'):
            seq = ""
            seq += line[1:].strip('translation="')
            while not line.endswith('"'):
                line = next(file_input).strip()
                seq += line
    return gene, seq