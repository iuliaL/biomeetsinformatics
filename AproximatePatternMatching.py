# Approximate Pattern Matching Problem:
# Find all approximate occurrences of a pattern in a string. See PatternMatchingGenome.py

# Input: Strings Pattern and Text along with an integer d which is the HammimgDistance between Pattern and Pattern'
# Output: All starting positions where Pattern appears as a substring of Genome with at most d mismatches. (A list of positions)

from HammingDistance import HammingDistance


def ApproximatePatternMatching(Pattern, Genome, d):
    positions = []
    k = len(Pattern)
    for index in range(len(Genome) - k + 1):
        curr_kmer = Genome[index : index + k]
        if HammingDistance(Pattern, curr_kmer) <= d:
            positions.append(index)

    return positions

# Test
# print(ApproximatePatternMatching('ATTCTGGA', 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT', 3)) # -> [6, 7, 26, 27]