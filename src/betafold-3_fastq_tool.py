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

