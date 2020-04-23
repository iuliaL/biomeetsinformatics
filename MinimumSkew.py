# Minimum Skew Problem: Find a position in a genome where the skew diagram attains a minimum.
# Input: A DNA string Genome.
# Output: All integer(s) i minimizing Skew[i] among all values of i (from 0 to len(Genome)).
from SkewArray import SkewArray


def MinimumSkew(Genome):
    skew = SkewArray(Genome)
    minimum = min(skew)
    min_positions = []
    for i in range(len(skew)):
        if skew[i] == minimum:
            min_positions.append(i)
    return min_positions

# => [3923620, 3923621, 3923622, 3923623]
