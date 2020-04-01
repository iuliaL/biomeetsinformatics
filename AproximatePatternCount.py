# Count how many times Pattern and slightly different Pattern (with d allowed differences) are found in Genome

from HammingDistance import HammingDistance

def ApproximatePatternCount(Pattern, Genome, d):
    k = len(Pattern)
    count = 0
    for index in range(len(Genome)):
        curr_kmer = Genome[index:index + k]
        if not len(curr_kmer) == k: # Hamming Distance assumes 2 equal length strings
            continue
        elif curr_kmer == Pattern or HammingDistance(Pattern, curr_kmer) <= d:
            count += 1 
    return count 


print(ApproximatePatternCount('AAAAA', 'AACAAGCATAAACATTAAAGAG', 1))  # => 4