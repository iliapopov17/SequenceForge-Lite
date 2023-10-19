complement_DNA = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
complement_RNA = {'A': 'U', 'G': 'C', 'U': 'A', 'C': 'G',
                'a': 'u', 'g': 'c', 'u': 'a', 'c': 'g',}


def transcribe(sequence: str) -> str:
    '''
    This function transcribes a DNA sequence to an RNA sequence

    Parameters:
    - sequence (str): DNA sequence

    Returns:
    - str: transcribed RNA sequence
    '''

    rna_sequence = ''
    for base in sequence:
        if base == 'T':
            rna_sequence += 'U'
        elif base == 't':
            rna_sequence += 'u'
        else:
            rna_sequence += base
    return rna_sequence


def reverse(sequence: str) -> str:
    '''
    This function reverses a DNA or RNA sequence

    Parameters:
    - sequence (str): DNA or RNA sequence

    Returns:
    - str: reversed sequence
    '''

    return sequence[::-1]


def complement(sequence: str) -> str:
    '''
    This function generates the complement of a DNA or RNA sequence

    Parameters:
    - sequence (str): DNA or RNA sequence

    Returns:
    - str: complemented sequence
    '''

    if ('T' in sequence) or ('t' in sequence):
        complement_sequence = ''.join([complement_DNA[nt] for nt in sequence])
    elif ('U' in sequence) or ('u' in sequence):
        complement_sequence = ''.join([complement_RNA[nt] for nt in sequence])

    return complement_sequence


def reverse_complement(sequence: str) -> str:
    '''
    This function generates the reverse complement of a DNA or RNA sequence

    Parameters:
    - sequence (str): DNA or RNA sequence

    Returns:
    - str: reverse complemented sequence
    '''

    return complement(reverse(sequence))


def run_dna_rna_tools(*args: str) -> str:
    '''
    Perform DNA or RNA operations based on the specified procedure

    Parameters:
    - *args (str):
        1. sequence(s) - one or multiple DNA or RNA sequences
        2. procedure - DNA or RNA processing procedure

    Returns:
    - str: result of the specified procedure on the sequences.
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