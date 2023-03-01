def rev_comp(sequence):
    complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    seq = sequence.upper()
    rcomp_seq = ''
    for base in seq[::-1]:
        rcomp_seq += complements[base]
    return rcomp_seq
