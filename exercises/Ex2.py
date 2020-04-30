# Compute the probability that ten randomly selected 15-mers from 10  600-nucleotide long strings 
# capture AT LEAST ONE implanted 15-mer.

"""Using one DNA string example: ttACCTtaac
(out of motifs with implanted ACGT)

The string has 10 characters and it is possible to extract 7 4-mers from it ( ttAC, tACC, ACCT, CCTt, CTta, Ttaa , taac).
Only one of these 4-mers are the implanted one, so the probability of randomly selecting it is 1/7.
The complementary probability 6/7 of missing the desired 4-mer.
Considering we have 5 strings, the probability of finding at least one implanted 4-mer is P = 1 - ((6/7)^5)"""

# Number of 15mers in 600 nucleotides
kmer_count = 600 - 15 + 1

# Probability that the kmer is the implanted one

pr_per_string = 1 / kmer_count 
pr_to_not_pick_correctly = 1 - pr_per_string

total_pr_to_not_pick_correctly = pr_to_not_pick_correctly ** 10
total_probability_to_pick_implanted = 1 - total_pr_to_not_pick_correctly
print(total_probability_to_pick_implanted)



