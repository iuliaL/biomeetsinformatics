# SkewArray's point is to calculate G-C occurences (so far) for each position of the Genome
# The final purpose is to find OriC for CIRCULAR baterial genome
# hipothesis: G - C
# positive on the forward half-strand (Ori -> term)
# negative on the reverse half-strand (term -> Ori)


import matplotlib.pyplot as plt


def calculateDiff(sy):
   return  -1 if sy == 'C' else 1 if sy == 'G' else 0

# Input:  A String Genome
# Output: The skew array of Genome as a list.

def SkewArray(Genome):
    skew = [0] # set zero position as zero on purpose
    for index in range(1, len(Genome) +  1):
        skew.append(skew[index - 1] + calculateDiff(Genome[index - 1]))
    return skew


if __name__ == "__main__":
    import sys
    if len(sys.argv[1:]) != 1:
        print("You must pass a genome string argument")
        exit()
    with open(sys.argv[1],'r') as file:
        genome = file.read()

        skew = SkewArray(genome)
        plt.plot(skew)
        plt.xlabel('Genome position')
        plt.ylabel('skew (difference G - C)')
        plt.show()



