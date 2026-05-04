import numpy as np

# WARNING: Unsorted
def find_sub(str, sub):
    """Returns all indices of sub in str"""
    start = 0
    while True:
        start = str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)


def site_instances(seq, site_seq):
    """Returns all indices of a restriction site sequence in seq. Only use with valid palindromic restriction site
    sequences"""
    # Concatenate forwards, and then sort
    instances = list(find_sub(seq, site_seq))
    instances.sort()
    return instances


def load_ref_by_id(fasta_iterator, id):
    """For reference genomes in the FNA format"""
    for fasta in fasta_iterator:
        if fasta.id == id:
            return fasta


def remove_outside(str, ranges):
    result = ''
    for range in ranges:
        result += str[range[0]:range[1]]
    return result