import os
from typing import Union, Tuple


from src import betafold_3_amino_analyzer
from src import betafold_3_dna_rna_tools
from src import betafold_3_fastq_tool


def run_amino_analyzer(sequence: str, procedure: str, *, weight_type: str = 'average', enzyme: str = 'trypsin'):
    '''
    This is the main function to run the amino-analyzer.py tool
    
    Parameters:
    - sequence (str): amino acid sequence in one-letter code
    - procedure (str): amino-analyzer.py tool has 5 functions at all:
        1. aa_weight - Calculates the amino acids weight in a protein sequence. Return float weight
            weight_type = 'average': default argument for 'aa_weight' function
            weight_type = 'monoisotopic' can be used as a second option
        2. count_hydroaffinity - Counts the quantity of hydrophobic and hydrophilic amino acids in a protein sequence. Return list in order: hydrophobic, hydrophilic
        3. peptide_cutter - Identifies cleavage sites in a given peptide sequence using a specified enzyme. Return list of cleavage sites
            enzyme = 'trypsin': default argument for 'peptide_cutter' function
            enzyme = 'chymotrypsin' can be used as a second option
        4. one_to_three_letter_code - Converts a protein sequence from one-letter amino acid code to three-letter code. Return string of amino acids in three-letter code
        5. sulphur_containing_aa_counter - Counts sulphur-containing amino acids in a protein sequence. Return quantaty of sulphur-containing amino acids

    Returns:
    - The result of the specified procedure
    '''

    procedures = ['aa_weight', 'count_hydroaffinity', 'peptide_cutter', 'one_to_three_letter_code', 'sulphur_containing_aa_counter']
    if procedure not in procedures:
        raise ValueError(f'Incorrect procedure. Acceptable procedures: {", ".join(procedures)}')

    if not betafold_3_amino_analyzer.is_aa(sequence):
        raise ValueError('Incorrect sequence. Only amino acids are allowed (V, I, L, E, Q, D, N, H, W, F, Y, R, K, S, T, M, A, G, P, C, v, i, l, e, q, d, n, h, w, f, y, r, k, s, t, m, a, g, p, c).')

    if procedure == 'aa_weight':
        result = betafold_3_amino_analyzer.aa_weight(sequence, weight_type)
    elif procedure == 'count_hydroaffinity':
        result = betafold_3_amino_analyzer.count_hydroaffinity(sequence)
    elif procedure == 'peptide_cutter':
        result = betafold_3_amino_analyzer.peptide_cutter(sequence, enzyme)
    elif procedure == 'one_to_three_letter_code':
        result = betafold_3_amino_analyzer.one_to_three_letter_code(sequence)
    elif procedure == 'sulphur_containing_aa_counter':
        result = betafold_3_amino_analyzer.sulphur_containing_aa_counter(sequence)
    return result


def run_dna_rna_tools(*args: str) -> str:
    '''
    This is the main function that performs DNA or RNA operations based on the specified procedure

    Parameters:
    - *args (str):
        1. sequence(s) - one or multiple DNA or RNA sequences
        2. procedure - DNA or RNA processing procedure

    Returns:
    - str: result of the specified procedure on the sequences
    '''

    if len(args) < 2:
        raise ValueError('Insufficient arguments. Indicate sequences and procedure.')

    sequences = args[:-1]
    procedure = args[-1]

    procedures = ['transcribe', 'reverse', 'complement', 'reverse_complement']
    if procedure not in procedures:
        raise ValueError(f'Incorrect procedure. Acceptable procedures: {", ".join(procedures)}')

    valid_bases = {'A', 'T', 'C', 'G', 'a', 't', 'c', 'g'}
    for sequence in sequences:
        if not set(sequence).issubset(valid_bases):
            raise ValueError('Incorrect sequence. Only nucleic acids are allowed (A, T, C, G, a, t, c, g).')

    results = []

    for sequence in sequences:
        if procedure == 'transcribe':
            result = betafold_3_dna_rna_tools.transcribe(sequence)
        elif procedure == 'reverse':
            result = betafold_3_dna_rna_tools.reverse(sequence)
        elif procedure == 'complement':
            result = betafold_3_dna_rna_tools.complement(sequence)
        elif procedure == 'reverse_complement':
            result = betafold_3_dna_rna_tools.reverse_complement(sequence)
        results.append(result)

    if len(results) == 1:
        return results[0]
    else:
        return results
    

def run_fastq_tool(
    input_path: str,
    gc_bounds: Union[int, float, Tuple[Union[int, float], Union[int, float]]]= (0, 100),
    length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
    quality_threshold: int = 0,
    output_filename: str = None) -> None:
    '''
    This is the main function that filters the DNA sequences in a fastq file based on specified criteria

    Parameters:
    - input_path (str): path to fastq file
    - gc_bounds (Union[int, float, Tuple[Union[int,float], Union[int,float]]): bottom and top bounds for GC content
    - length_bounds (Union[int, Tuple[int, int]]): bottom and top bounds for sequence length
    - quality_threshold (int): acceptable quality
    - output_filename (str): name of filtered fastq file

    Saves the filtered sequences to a fastq file named input_path/output_filename
    '''
    seqs = betafold_3_fastq_tool.read_fastq(input_path)

    if isinstance(gc_bounds, int) or isinstance(gc_bounds, float):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int) or isinstance(length_bounds, float):
        length_bounds = (0, length_bounds)

    filtered_seqs = {}
    
    for name, (sequence, quality) in seqs.items():        
        if (betafold_3_fastq_tool.gc_check(sequence, gc_bounds) 
            and betafold_3_fastq_tool.length_check(sequence, length_bounds) 
            and betafold_3_fastq_tool.quality_check(quality, quality_threshold)):
            filtered_seqs[name] = (sequence, quality)
    
    if output_filename is None:
        output_filename = os.path.basename(input_path)
        
    else:
        output_filename += '.fastq'
        
    betafold_3_fastq_tool.write_fastq(output_filename, filtered_seqs)