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


def rc_naive(pattern, sequence):
    """ This function performs a naive exact match alignment of pattern and the reverse
        complement of pattern onto sequence

    Args:
        pattern (str): A pattern to attempt to align to sequence
        sequence (str): The sequence to attempt alignment to.

    Returns:
        list [ ints ]: The list of indices where pattern can be exactly aligned to sequence
    """
    matches = []
    pattern_rc = rev_comp(pattern)
    for sidx in range(len(sequence) - len(pattern) + 1):
        matched = True
        for pidx in range(len(pattern)):
            if sequence[sidx + pidx] != pattern[pidx]:
                matched = False
                break
        if matched:
            matches.append(sidx)
    for sidx in range(len(sequence) - len(pattern_rc) + 1):
        matched = True
        for pidx in range(len(pattern_rc)):
            if sequence[sidx + pidx] != pattern_rc[pidx]:
                matched = False
                break
        if matched and sidx not in matches:
            matches.append(sidx)
    return matches


def seq_from_fasta(filename):
    """ Strips the initial line from a fasta file containing a single sequence
        for a single organism. Reads the remaining lines and concatenates them
        onto an str.

    Args:
        filename (str): The name of the file to read from in fasta format

    Returns:
        str: The genome sequence contained in the file
    """
    sequence = ''
    try:
        file = open(filename)
    except FileNotFoundError:
        print(f"{filename} not found!")
        return "Error, file not found"
    for line in file:
        line = line.rstrip()
        if not line.startswith('>'):
            sequence += line
    return sequence


def naive_2mm(pattern, sequence):
    """ This function performs a naive exact match alignment of pattern onto sequence
        allowing for 2 mismatches per occurence

    Args:
        pattern (str): A pattern to attempt to align to sequence
        sequence (str): The sequence to attempt alignment to.

    Returns:
        list [ ints ]: The list of indices where pattern can be aligned to sequence
    """
    matches = []
    for sidx in range(len(sequence) - len(pattern) + 1):
        mm = 0
        matched = True
        for pidx in range(len(pattern)):
            if sequence[sidx + pidx] != pattern[pidx]:
                mm += 1
                if mm > 2:
                    matched = False
                    break
        if matched:
            matches.append(sidx)
    return matches


def reads_from_fastq(filename):
    """ Strips the initial line from a fastq file containing a single sequence
        for a single organism. Reads the remaining lines and concatenates them
        onto an str.

    Args:
        filename (str): The name of the file to read from in fastq format

    Returns:
        [str] [str] The sequences and the quality values encoded with phred33
    """
    seqs = []
    quals = []
    try:
        file = open(filename)
    except FileNotFoundError:
        print(f"{filename} not found!")
        return "Error, file not found"
    while True:
        file.readline()
        seq = file.readline().rstrip()
        file.readline()
        qual = file.readline().rstrip()
        if len(seq) == 0:
            break
        seqs.append(seq)
        quals.append(qual)
    return seqs, quals
