def rev_comp(sequence):
    """ rev_comp returns the reverse complement to a sequence of nucleotides. The complement
        of a base is stored in a dictionary look up table, with N being its own complement.
        The reverse complement of a sequence is the sequence in reverse order with all bases
        replaced by their complements. Useful in genetic sequencing.

    Args:
        sequence (str): The sequence of nucleotides to compute the reverse complement

    Returns:
        str: The reverse complement of sequence
    """
    complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    seq = sequence.upper()
    rcomp_seq = ''
    for base in seq[::-1]:
        rcomp_seq += complements[base]
    return rcomp_seq


def naive(pattern, sequence):
    """ This function performs a naive exact match alignment of pattern onto sequence

    Args:
        pattern (str): A pattern to attempt to align to sequence
        sequence (str): The sequence to attempt alignment to.

    Returns:
        list [ ints ]: The list of indices where pattern can be exactly aligned to sequence
    """
    matches = []
    for sidx in range(len(sequence) - len(pattern) + 1):
        matched = True
        for pidx in range(len(pattern)):
            if sequence[sidx + pidx] != pattern[pidx]:
                matched = False
                break
        if matched:
            matches.append(sidx)
    return matches
