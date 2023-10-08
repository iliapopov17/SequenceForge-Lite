def gc_check(sequence: str, gc_bounds: tuple) -> bool:
    """
    This function checks if the GC content of a DNA sequence in fastq file lies within specified bounds

    Parameters:
    - seq (str): DNA sequence
    - gc_bounds (tuple): tuple with bottom and top bounds for GC content

    Returns:
    - bool: True if GC content is within bounds, False if not
    """
    
    gc = 0

    for base in sequence:
        if base == "G" or base == "C":
            gc += 1
            
    bottom_threshold, top_threshold = gc_bounds

    return bottom_threshold <= gc <= top_threshold

def length_check(sequence: str, length_bounds: tuple) -> bool:
    """
    This function checks if the length of a DNA sequence in fastq file lies within specified bounds

    Parameters:
    - seq (str): DNA sequence
    - length_bounds (tuple): tuple with bottom and top bounds for sequence length

    Returns:
    - bool: True if the length is within bounds, False if not
    """
    
    bottom_threshold, top_threshold = length_bounds

    return bottom_threshold <= len(sequence) <= top_threshold


def quality_check(sequence: str, quality_threshold: int) -> bool:
    """
    This function checks if the quality score of a DNA sequence in fastq file lies above a specified bound.

    Parameters:
    - seq (str): DNA sequence
    - quality_threshold (int): acceptable quality

    Returns:
    - bool: True if the average quality score is above the threshold, False if not
    """
    
    quality = 0

    for char in sequence:
        quality += ord(char) - 33

    return quality / len(sequence) >= quality_threshold