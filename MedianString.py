from HammingDistance import HammingDistance
from Combinations import Combinations

# this function might not be right
def MedianString(k, *Dna):
    median = None
    for comb in Combinations(k):
        for string in Dna:
            distance = float('inf')
            for i in range(len(string) - k + 1):
                curr_kmer = string[i: i + k]
                curr_distance = HammingDistance(curr_kmer, comb)
                if distance > curr_distance:
                    median = comb
                    distance = curr_distance
    return median


