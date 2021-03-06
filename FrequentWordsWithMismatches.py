# Find the most frequent kmers with at most d mismatches

from AproximatePatternCount import ApproximatePatternCount
from Neighbors import Neighbors
from Combinations import Combinations

import itertools
import time
from collections import defaultdict


def FrequentWordsWithMismatches(Genome, k, d):
    start = time.process_time()

    '''
    1. Create a k-mer list with ['A','C','G','T'] (4**k k-mers in the list)
    2. Calculate approximate pattern count for each k-mer (total 4**k)
    3. Get k-mer(s) with maximum approximate pattern count.
    '''
    aprox_frq_words = []
    frequencies = {}
    # all existent combinations of k-length kmers 4**k
    combinations = Combinations(k)
    for kmer in combinations:
        count = ApproximatePatternCount(kmer, Genome, d)
        frequencies[kmer] = count
    for kmer in frequencies:
        if frequencies[kmer] == max(frequencies.values()):
            aprox_frq_words.append(kmer)
    end = time.process_time()
    print("Time:", end - start)
    return aprox_frq_words

def FasterFrequentWordsWithMismatches(Genome, k, d):
    start = time.process_time()
    aprox_frq_words = []
    frequencies = defaultdict(lambda: 0)
    # all existent kmers with d mismatches of current kmer in genome
    for index in range(len(Genome) - k + 1):
        curr_kmer = Genome[index : index + k]
        neighbors = Neighbors(curr_kmer, d)
        for kmer in neighbors:
            frequencies[kmer] += 1 

    for kmer in frequencies:
        if frequencies[kmer] == max(frequencies.values()):
            aprox_frq_words.append(kmer)
    end = time.process_time()
    print("Time:", end - start)
    return aprox_frq_words




# print(FrequentWordsWithMismatches("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4,1))
# print(FasterFrequentWordsWithMismatches("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4,1))



