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
    # skewAsDict =  { i : skew[i] for i in range(0, len(skew) ) }
    # return skewAsDict
    # print(skew[:])
    return skew


with open('E_coli.txt') as file:
    e_coli = file.read()

    skew_E_Coli = SkewArray(e_coli)

    # plt.plot(*zip(*skew_E_Coli.items()))
    plt.plot(skew_E_Coli,  marker='o')
    plt.show()



print(SkewArray("CATGGGCATCGGCCATACGCC"))