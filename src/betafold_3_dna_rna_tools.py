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