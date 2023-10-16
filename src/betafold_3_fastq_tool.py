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
    
    gc_percent = gc / len(sequence) * 100
            
    bottom_threshold, top_threshold = gc_bounds

    return bottom_threshold <= gc_percent <= top_threshold


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


def run_fastq_tool(
    seqs: dict,
    gc_bounds: tuple = (0, 100),
    length_bounds: tuple = (0, 2**32),
    quality_threshold: int = 0,
) -> dict:
    """
    This is the main function that filters the DNA sequences in a fastq file based on specified criteria

    Parameters:
    - seqs (dict): dictionary where keys are sequence identifiers and values are tuples containing sequence and quality scores
    - gc_bounds (tuple): tuple with bottom and top bounds for GC content
    - length_bounds (tuple): tuple with bottom and top bounds for sequence length
    - quality_threshold (int): acceptable quality

    Returns:
    - dict: filtered dictionary that contains only the sequences that fully match the specified criteria
    """

    # ------------------------------------------------------
    # Check if every parameter is setted up correctly
    if type(gc_bounds) == int:
        gc_bounds = (0, gc_bounds)

    elif len(gc_bounds) == 1:
        gc_bounds = (0, gc_bounds[0])

    elif len(gc_bounds) > 2:
        raise ValueError("Invalid gc_bound value!")

    elif type(gc_bounds) != tuple:
        raise ValueError("Input tuple value!")

    if type(length_bounds) == int:
        length_bounds = (0, length_bounds)

    elif len(length_bounds) == 1:
        length_bounds = (0, length_bounds[0])

    elif len(length_bounds) > 2:
        raise ValueError("Invalid length_bounds value!")

    elif type(length_bounds) != tuple:
        raise ValueError("Input tuple value!")

    if len(str(quality_threshold)) > 2:
        raise ValueError("Invalid quality_threshold value!")
    
    elif type(quality_threshold) != int:
        raise ValueError("Input integer value!")
    # ------------------------------------------------------

    filtered_seqs = {}
    
    for name, (sequence, quality) in seqs.items():        
        if ((gc_check(sequence, gc_bounds))
            and (length_check(sequence, length_bounds))
            and (quality_check(quality, quality_threshold))):
            filtered_seqs[name] = (sequence, quality)
    
    return filtered_seqs