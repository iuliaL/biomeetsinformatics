# Find the most frequent k-mers (with mismatches and reverse complements) in a string.

import time
from collections import defaultdict

from Neighbors import Neighbors
from ReverseComplement import ReverseComplement

def FrequentWordsWithMismatchesAndReverseComplements(Genome, k, d):
    start = time.process_time()
    aprox_frq_words = []
    frequencies = defaultdict(lambda: 0)
    # all existent kmers with d mismatches of current kmer in genome
    for index in range(len(Genome) - k + 1):
        curr_kmer = Genome[index : index + k]
        neighbors = Neighbors(curr_kmer, d)
        for kmer in neighbors:
            frequencies[kmer] += 1
            frequencies[ReverseComplement(kmer)] += 1

    for kmer in frequencies:
        if frequencies[kmer] == max(frequencies.values()):
            aprox_frq_words.append(kmer)
    end = time.process_time()
    print("Time:", end - start)
    return aprox_frq_words


# print(FrequentWordsWithMismatchesAndReverseComplements("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1))



