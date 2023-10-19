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

    if not is_aa(sequence):
        raise ValueError('Incorrect sequence. Only amino acids are allowed (V, I, L, E, Q, D, N, H, W, F, Y, R, K, S, T, M, A, G, P, C, v, i, l, e, q, d, n, h, w, f, y, r, k, s, t, m, a, g, p, c).')

    if procedure == 'aa_weight':
        result = aa_weight(sequence, weight_type)
    elif procedure == 'count_hydroaffinity':
        result = count_hydroaffinity(sequence)
    elif procedure == 'peptide_cutter':
        result = peptide_cutter(sequence, enzyme)
    elif procedure == 'one_to_three_letter_code':
        result = one_to_three_letter_code(sequence)
    elif procedure == 'sulphur_containing_aa_counter':
        result = sulphur_containing_aa_counter(sequence)
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
            result = transcribe(sequence)
        elif procedure == 'reverse':
            result = reverse(sequence)
        elif procedure == 'complement':
            result = complement(sequence)
        elif procedure == 'reverse_complement':
            result = reverse_complement(sequence)
        results.append(result)

    if len(results) == 1:
        return results[0]
    else:
        return results
    

def run_fastq_tool(
    seqs: dict,
    gc_bounds: tuple = (0, 100),
    length_bounds: tuple = (0, 2**32),
    quality_threshold: int = 0,
) -> dict:
    '''
    This is the main function that filters the DNA sequences in a fastq file based on specified criteria

    Parameters:
    - seqs (dict): dictionary where keys are sequence identifiers and values are tuples containing sequence and quality scores
    - gc_bounds (tuple): tuple with bottom and top bounds for GC content
    - length_bounds (tuple): tuple with bottom and top bounds for sequence length
    - quality_threshold (int): acceptable quality

    Returns:
    - dict: filtered dictionary that contains only the sequences that fully match the specified criteria
    '''

    # ------------------------------------------------------
    # Check if every parameter is setted up correctly
    if type(gc_bounds) == int:
        gc_bounds = (0, gc_bounds)

    elif len(gc_bounds) == 1:
        gc_bounds = (0, gc_bounds[0])

    elif len(gc_bounds) > 2:
        raise ValueError('Invalid gc_bound value!')

    elif type(gc_bounds) != tuple:
        raise ValueError('Input tuple value!')

    if type(length_bounds) == int:
        length_bounds = (0, length_bounds)

    elif len(length_bounds) == 1:
        length_bounds = (0, length_bounds[0])

    elif len(length_bounds) > 2:
        raise ValueError('Invalid length_bounds value!')

    elif type(length_bounds) != tuple:
        raise ValueError('Input tuple value!')

    if len(str(quality_threshold)) > 2:
        raise ValueError('Invalid quality_threshold value!')
    
    elif type(quality_threshold) != int:
        raise ValueError('Input integer value!')
    # ------------------------------------------------------

    filtered_seqs = {}
    
    for name, (sequence, quality) in seqs.items():        
        if (gc_check(sequence, gc_bounds) 
            and length_check(sequence, length_bounds) 
            and quality_check(quality, quality_threshold)):
            filtered_seqs[name] = (sequence, quality)
    
    return filtered_seqs