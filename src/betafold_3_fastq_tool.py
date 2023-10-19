from typing import Dict, Tuple


def read_fastq(input_path: str) -> Dict[str, Tuple[str, str]]:
    '''
    This function reads fastq file and makes from dict of sequences

    Parameters:
    - input_path (str): path to a fastq file

    Returns:
    - Dict: dictionary:
        1. keys - seq names
        2. values - tuple of sequences and their's quality
    '''

    seqs = {}
    with open(input_path, 'r') as file:
        for line in file:
            if line.startswith("@"):
                counter = 0
                name = line[1:]
            if counter == 1:
                seq = line
            if counter == 3:
                quality = line
                seqs[name] = (seq, quality)
            counter += 1
    return seqs


def write_fastq(output_filename: str, seqs: Dict[str, Tuple[str, str]]):
    '''
    This function writes a dictionary:
        1. keys - sequence names
        2. values - a tuple of a sequence its quality
    to a file in fastq format

    Parameters:
    - output_filename (str): name of an output fastq file
    - seqs (Dict[str, Tuple[str, str]]): dictionary of sequences
    '''

    with open(output_filename, 'w') as file:
        for name, seq in seqs.items():
            file.write(f'@{name}')
            file.write(f'{seq[0]}')
            file.write(f'+{name}')
            file.write(f'{seq[1]}')


def gc_check(sequence: str, gc_bounds: tuple) -> bool:
    '''
    This function checks if the GC content of a DNA sequence in fastq file lies within specified bounds

    Parameters:
    - seq (str): DNA sequence
    - gc_bounds (tuple): tuple with bottom and top bounds for GC content

    Returns:
    - bool: True if GC content is within bounds, False if not
    '''
    
    gc = 0
    sequence = sequence.upper()

    for base in sequence:
        if base == 'G' or base == 'C':
            gc += 1
    
    gc_percent = gc / len(sequence) * 100
            
    bottom_threshold, top_threshold = gc_bounds

    return bottom_threshold <= gc_percent <= top_threshold


def length_check(sequence: str, length_bounds: tuple) -> bool:
    '''
    This function checks if the length of a DNA sequence in fastq file lies within specified bounds

    Parameters:
    - seq (str): DNA sequence
    - length_bounds (tuple): tuple with bottom and top bounds for sequence length

    Returns:
    - bool: True if the length is within bounds, False if not
    '''
    
    bottom_threshold, top_threshold = length_bounds

    return bottom_threshold <= len(sequence) <= top_threshold


def quality_check(sequence: str, quality_threshold: int) -> bool:
    '''
    This function checks if the quality score of a DNA sequence in fastq file lies above a specified bound.

    Parameters:
    - seq (str): DNA sequence
    - quality_threshold (int): acceptable quality

    Returns:
    - bool: True if the average quality score is above the threshold, False if not
    '''
    
    quality = 0

    for char in sequence:
        quality += ord(char) - 33

    return quality / len(sequence) >= quality_threshold