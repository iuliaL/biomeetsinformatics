# Pattern Matching Problem: Find all occurrences of a pattern in a string.
#     Input: Strings Pattern and Genome.
#     Output: All starting positions in Genome where Pattern appears as a substring.


def PatternMatching(Pattern, Genome):
    positions = []
    k = len(Pattern)
    for index in range(len(Genome) - k + 1):
        currentKmer = Genome[index:index + k]
        if currentKmer == Pattern:
            positions.append(index)
    return positions

pattern = "CTTGATCAT"
genomeFile = open(r"vibrio_cholerae_genome.txt")
genome = genomeFile.read()
genomeFile.close()


print("The positions where I found the pattern CTTGATCAT in V.Cholerae genome are: " +
', '.join(
    map(str, PatternMatching(pattern, genome)))
    )
"""
The positions where I found the pattern CTTGATCAT in V.Cholerae genome are:
60039, 98409, 129189, 152283, 152354, 152411, 163207, 197028, 200160, 357976, 376771, 392723, 532935, 600085, 622755, 1065555

152283, 152354, 152411 form a "clump" a.k.a appear close to each other in a small region of the genome.

"""
