# What is the expected number of occurrences of a 9-mer in 500 random DNA strings, each of length 1000? Assume that the sequences are formed by selecting each nucleotide (A, C, G, T) with the same probability (0.25).

# In order to answer first we have to think of the following:

# 1.How many 9-mers exists? (this will not enter in the calculus)
4 ** 9 -> 262144 # ( 4 nucleotides, 4**k )

# 2.What's the probability of, in a set of the 9-mers of Hint 1, you pick a 9-mer that you want in the first pick?
1 / 262144

# 3.How many 9-mers exists in a string with 1000 nucleotides?
1000 - 9 + 1 = 992

# 4.How many 9-mers exists in a set of 500 strings of length 1000?
992 * 500 = 496000

# 5.What is the expected number of occurrences of the 9-mer that you choose in Hint 2 occur in the set of Hint 4?
:)

262144 ..... 1
496000 ......x

496000 /  262144
