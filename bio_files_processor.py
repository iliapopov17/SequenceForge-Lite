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