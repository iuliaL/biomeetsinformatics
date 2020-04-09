from HammingDistance import HammingDistance
from random import randint


def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for key in "ACGT":
        count[key] = [0] * k  # create a matrix of 4 rows by k length filled with zeros
    for motif in Motifs:
        for index in range(k):
            count[motif[index]][index] += 1
    return count

# print("Count Matrix", Count(['AACGTA','CCCGTT', 'CACCTT', 'GGATTA','TTCCGG']))


def Profile(Motifs):  # this is just like a percentage
    count = Count(Motifs)
    for nucleotide_key in count:
        count[nucleotide_key] = [x / len(Motifs) for x in count[nucleotide_key]]
    return count


def CountWithPseudocounts(Motifs):
    count = {}
    k = len(Motifs[0])
    pseudocount = 1  # this adds 1 for each nucleotide count in order to avoid computing zero probabilities later
    for key in "ACGT":
        count[key] = [pseudocount] * k  # create a matrix of 4 rows by k length filled with 1
    for motif in Motifs:
        for index in range(k):
            count[motif[index]][index] += 1
    return count


def ProfileWithPseudocounts(Motifs):  # this is just like a percentage
    count = CountWithPseudocounts(Motifs)
    pseudocounts = 1 * 4  # (4 because it's 1 for each nucleotide)
    for nucleotide_key in count:
        count[nucleotide_key] = [x / (len(Motifs) + pseudocounts) for x in count[nucleotide_key]]
    return count


def Consensus(Motifs):
    consensus = ''
    count = CountWithPseudocounts(Motifs)
    k = len(Motifs[0])
    for index in range(k):
        max_so_far = 0
        most_freq_nucleotide = None
        for nucleotide, val in count.items():
            if val[index] > max_so_far:
                max_so_far = val[index]
                most_freq_nucleotide = nucleotide
        consensus += most_freq_nucleotide
    return consensus

# print("Consensus", Consensus(['AACGTA','CCCGTT', 'CACCTT', 'GGATTA','TTCCGG'])) # => CACCTA


def Score(Motifs):  # column by column distance to consensus
    consensus = Consensus(Motifs)
    count = CountWithPseudocounts(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    score = 0
    for index in range(k):
        subtract = 0
        nucleotide = consensus[index]
        for key, val in count.items():
            if key == nucleotide:
                subtract += val[index]
        score += t - subtract
    return score


def Score_(Motifs):  # row by row distance to consensus
    score = 0
    consensus = Consensus(Motifs)
    for motif in Motifs:
        score += HammingDistance(motif, consensus)
    return score

# print(Score(['AACGTA','CCCGTT', 'CACCTT', 'GGATTA','TTCCGG']))

# Probability


def Pr(Pattern, Profile):
    pr = 1
    for index in range(len(Pattern)):
        nucleotide = Pattern[index]
        pr *= Profile[nucleotide][index]
    return pr


# profile = {
#     'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
#     'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
#     'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
#     'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
# }

# print("Probability", Pr("TCGTGGATTTCC", profile))


def ProfileMostProbableKmer(text, k, profile):  # but here i could get more than 1, i ignore ties
    most_probable = None
    initial_probability = -1  # impossible one
    for i in range((len(text) - k + 1)):
        kmer = text[i:i + k]
        probability = Pr(kmer, profile)
        if probability > initial_probability:
            initial_probability = probability
            most_probable = kmer
    return most_probable


# profile__ = {
#     'A': [0.2, 0.2, 0.3, 0.2, 0.3],
#     'C': [0.4, 0.3, 0.1, 0.5, 0.1],
#     'G': [0.3, 0.3, 0.5, 0.2, 0.4],
#     'T': [0.1, 0.2, 0.1, 0.1, 0.2]
# }

# print("Most probable by profile", ProfileMostProbableKmer("ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT", 5, profile__))    # => CCGAG


# http://www.mrgraeme.co.uk/greedy-motif-search/
# Amazing explanation

def GreedyMotifSearch(Dna, k, t):  # Dna is a list of t strings (dont know why i have t since it's len(Dna))
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    dna_string_length = len(Dna[0])
    for i in range(dna_string_length - k + 1):
        Motifs = []  # initial Motifs
        kmer_in1st_string = Dna[0][i: i + k]
        Motifs.append(kmer_in1st_string)
        for index in range(1, t):
            profile_matrix = ProfileWithPseudocounts(Motifs)
            most_similar = ProfileMostProbableKmer(Dna[index], k, profile_matrix)
            Motifs.append(most_similar)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


# Dna = [
#     "TTACCTTAAC",
#     "GATGTCTGTC",
#     "ACGGCGTTAG",
#     "CCCTAACGAG",
#     "CGTCAGAGGT"]

# print("Greedy motif search", GreedyMotifSearch(Dna, 4, 5))  # => ['ACCT', 'ATGT', 'ACGG', 'ACGA', 'AGGT']

# Dna_ = [
#     "GGCGTTCAGGCA",
#     "AAGAATCAGTCA",
#     "CAAGGAGTTCGC",
#     "CACGTCAATCAC",
#     "CAATAATATTCG"
# ]
# print("Greedy motif search", GreedyMotifSearch(Dna_, 3, 5))  # => ['CAG', 'CAG', 'CAA', 'CAA', 'CAA']


def Motifs(Profile, k, Dna):
    result = []
    for string in Dna:
        most_probable = ProfileMostProbableKmer(string, k, Profile)
        result.append(most_probable)
    return result


def RandomMotifs(Dna, k, t):
    random_motifs = []
    for s in Dna:
        start_pos = randint(0, len(Dna) - 1)
        random_motifs.append(s[start_pos: start_pos + k])
    return random_motifs


def RandomizedMotifSearch(Dna, k, t):
    BestMotifs = RandomMotifs(Dna, k, t)
    while True:
        profile = ProfileWithPseudocounts(BestMotifs)
        new_motifs = Motifs(profile, k, Dna)
        if Score(new_motifs) < Score(BestMotifs):
            BestMotifs = new_motifs
        else:
            return BestMotifs


Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",
"CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",
"ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",
"GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",
"GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",
"CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",
"GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",
"GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",
"GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",
"TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

# set t equal to the number of strings in Dna, k equal to 15, and N equal to 100.
t = len(Dna)
k = 15
N = 100

# Call RandomizedMotifSearch(Dna, k, t) N times, storing the best-scoring set of motifs
# resulting from this algorithm in a variable called BestMotifs
i = 0
BestMotifs = RandomizedMotifSearch(Dna, k, t)
while i < N:
    motifs = RandomizedMotifSearch(Dna, k, t)
    if Score(BestMotifs) > Score(motifs):
        BestMotifs = motifs
    i+=1

print(BestMotifs)
print(Score(BestMotifs))