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

pattern1 = "CTTGATCAT"
pattern2 = "ATGATCAAG" # reverse complement of pattern1
genomeFile = open(r"vibrio_cholerae_genome.txt")
genome = genomeFile.read()
genomeFile.close()


print("The positions where I found the pattern CTTGATCAT in V.Cholerae genome are: " +
', '.join(
    map(str, PatternMatching(pattern1, genome)))
    )

print("The positions where I found the pattern ATGATCAAG in V.Cholerae genome are: " +
', '.join(
    map(str, PatternMatching(pattern2, genome)))
    )
"""
The positions where I found the pattern CTTGATCAT in V.Cholerae genome are:
60039, 98409, 129189, 152283, 152354, 152411, 163207, 197028, 200160, 357976, 376771, 392723, 532935, 600085, 622755, 1065555

The positions where I found the pattern ATGATCAAG in V.Cholerae genome are: 
116556, 149355, 151913, 152013, 152394, 186189, 194276, 200076, 224527, 307692, 479770, 610980, 653338, 679985, 768828, 878903, 985368

Conclusion:
152283, 152354, 152411 for CTTGATCAT form a "clump" a.k.a appear close to each other in a small region of the genome.
idem for 151913, 152013, 152394 for ATGATCAAG

CTTGATCAT
||||||||
GAACTAGTA

Facts:
    ->  An Ori usually has ~500 nucleotides. For v.cholerae the number is 540
    ->  9-mers have been experimentally proven to be DnaA boxes
"""