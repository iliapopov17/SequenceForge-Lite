import os
from typing import List, Tuple

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


def select_gene_from_gbk(input_gbk: str, target_gene: str, n_before: int, n_after: int) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    '''
    This function selects genes around gene of interest in gbk file

    Parameters:
    - input_gbk (str): path to the gbk file
    - target_gene (str): name of the target gene
    - n_before (int): number of genes we want to get before target
    - n_after (int): number of genes we want to get after target

    Returns:
    - Tuple of two lists: genes before and after the target, consisting of tuples with the name of the gene and its amino acid sequence
    '''

    with open(input_gbk) as file_input:
        genes_before = [None] * n_before
        genes_after = [None] * n_after
        gene, seq = get_gene_sequence(file_input)
        while gene:
            if target_gene in gene:
                break
            genes_before.append((gene + f'_before_{target_gene}', seq))
            del genes_before[0]
            gene, seq = get_gene_sequence(file_input)
        while gene:
            genes_after.append((gene + f'_after_{target_gene}', seq))
            del genes_after[0]
            if None not in genes_after:
                return [g for g in genes_before if g], [g for g in genes_after if g]
    return [g for g in genes_before if g], [g for g in genes_after if g]


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: List[str], n_before: int, n_after: int, output_fasta: str = None) -> None:
    '''
    This function extracts gene sequences from the GenBank (gbk) file and creates a FASTA file with the specified neighbouring genes

    Parameters:
    - input_gbk (str): path to the input GenBank file
    - genes (List[str]): list of genes to extract
    - n_before (int): number of genes before the specified genes to include
    - n_after (int): number of genes after the specified genes to include
    - output_fasta (str): path to the output FASTA file

    Returns:
    - FASTA file
    '''

    if output_fasta is None:
        output_fasta = os.path.splitext(os.path.basename(input_gbk))[0] + '_target_proteins.fasta'
    with open(output_fasta, 'w') as file:
        for target in genes:
            genes_before, genes_after = select_gene_from_gbk(input_gbk, target, n_before, n_after)
            for gene, seq in genes_before:
                file.write(f'>{gene}\n{seq}\n')
            for gene, seq in genes_after:
                file.write(f'>{gene}\n{seq}\n')