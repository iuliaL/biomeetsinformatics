from itertools import product

# output a list of all possible kmers of k length

def Combinations(k, of='ACGT'):
    return ["".join(comb) for comb in list(product(of, repeat=k))]
