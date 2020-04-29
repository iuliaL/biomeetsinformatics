# Count how many times Pattern and slightly different Pattern (with d allowed differences) are found in Genome

from HammingDistance import HammingDistance

def ApproximatePatternCount(Pattern, Genome, d):
    k = len(Pattern)
    count = 0
    for index in range(len(Genome) - k + 1):
        curr_kmer = Genome[index:index + k]
        if HammingDistance(Pattern, curr_kmer) <= d:
            count += 1 
    return count 


# print(ApproximatePatternCount('ACAA', 'AACAAGCTGATAAACATTTAAAGAG', 1))  # => 4
