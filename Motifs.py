from HammingDistance import HammingDistance
from random import randint, uniform, random
from FrequentWords import FrequencyMap
from Neighbors import Neighbors
from AproximatePatternCount import ApproximatePatternCount

# This function should return a list of strings.
def MotifEnumeration(k, d, *dna):
    '''
        Patterns ← an empty set
        for each k-mer Pattern in the first string in Dna
            for each k-mer Pattern’ differing from Pattern by at most d mismatches
                if Pattern' appears in each string from Dna with at most d mismatches
                    add Pattern' to Patterns
        remove duplicates from Patterns
        return Patterns
    '''
    motifs = set()
    for i in range(len(dna[0]) - k + 1):
        kmer = dna[0][i: i + k]
        neighbors = Neighbors(kmer, d)
        for neighbor in neighbors:
            present = True
            index = 0
            while present and index < len(dna):
                string = dna[index]
                count = ApproximatePatternCount(neighbor, string, d)
                if count == 0:
                    present = False
                index +=  1
            if present:
                motifs.add(neighbor)
    return list(motifs)

# ??
# def FasterMotifEnumeration(k, d, *dna):
#     '''
#         for each k-mer Pattern in the strings of Dna
#         compute neighbors
#         and check intersections
#         return Patterns
#     '''
#     motif_sets = [set()] * len(dna)
#     for index, string in enumerate(dna):
#         for i in range(len(string) - k + 1):
#             kmer = string[i: i + k]
#             neighbors = Neighbors(kmer, d)
#             motif_sets[index].update(neighbors)
#     motif_sets =motif_sets
#     motifs =  set.intersection(*motif_sets)
#     return list(motifs)
# ??


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

# if there are multiple Profile-most probable k-mers in Text, then we select the first such k-mer occurring in Text.
def ProfileMostProbableKmer(text, k, profile):  # but here we could get more than 1, we ignore ties
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

# print("Most probable 5mer in sequence given a profile", ProfileMostProbableKmer("ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT", 5, profile__))
# => CCGAG


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


# Find the most probable kmer(motif) in each DNA string given the Profile
# Output the Motifs
def Motifs(Profile, k, Dna):
    motifs = []
    for string in Dna:
        most_probable = ProfileMostProbableKmer(string, k, Profile)
        motifs.append(most_probable)
    return motifs


def RandomMotifs(Dna, k, t):  # just pick 1 random k mer from each DNA string
    random_motifs = []
    for s in Dna:
        start_pos = randint(0, len(s) - k)
        random_motifs.append(s[start_pos: start_pos + k])
    return random_motifs

def RandomizedMotifSearch(Dna, k, t):
    BestMotifs = RandomMotifs(Dna, k, t)
    while True:
        # create a profile for the random motifs I start with and then use it to find better motifs and so on
        profile = ProfileWithPseudocounts(BestMotifs)
        # which runs a ProfileMostProbableKmer on EACH of the Dna strings, thus can change all motifs found before, good/bad
        new_motifs = Motifs(profile, k, Dna)
        if Score(new_motifs) < Score(BestMotifs):
            BestMotifs = new_motifs
        else:
            return BestMotifs


Dna = [
    "TTACCTTAAC",
    "GATGTCTGTC",
    "ACGGCGTTAG",
    "CCCTAACGAG",
    "CGTCAGAGGT"]

t = len(Dna)
k = 4
N = 100

times = []


def RunNTimesRandomizedMotifSearch(N):
    i = 0
    BestMotifs = RandomizedMotifSearch(Dna, k, t)
    while i < N:
        motifs = RandomizedMotifSearch(Dna, k, t)
        if Score(BestMotifs) > Score(motifs):
            BestMotifs = motifs
        i += 1
    consensus = Consensus(BestMotifs)
    times.append(consensus)


def Frequency(patterns):
    frequency = {}
    for Pattern in patterns:
        if Pattern not in frequency:
            frequency[Pattern] = 1
        else:
            frequency[Pattern] += 1
    return frequency


# i = 0
# while i < N:
#     RunNTimesRandomizedMotifSearch(N)
#     i += 1

# implanted motif (consensus) was ACGT indeed
# => {'ACGT': 83, 'TCAG': 3, 'GTTA': 2, 'CGTC': 3, 'GACG': 1, 'CCTT': 5, 'CTTA': 1, 'GCCT': 1, 'TTAG': 1}


# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities

def Normalize(Probabilities):
    total = sum(Probabilities.values())
    new_probabilities = {}
    for k in Probabilities:
        new_probabilities[k] = Probabilities[k] / float(total)
    return new_probabilities


# print(Normalize({'AC': 0.1, 'CG': 0.1, 'TG': 0.1, 'TT': 0.1}))

def WeightedDie(Probabilities): # returns one kmer from the probabilities dict
    n = uniform(0, 1)
    for p in Probabilities:
        n -= Probabilities[p]
        if n <= 0:
            return p

# print(WeightedDie({'AA': 0.3, 'AC': 0.2, 'TT': 0.45, 'CG': 0.05}))


def ProfileRandomlyGeneratedString(Text, profile, k):
    # Differs from ProfileMostProbableKmer because it is not necessarily picking the most probable kmer
    # (it has a degree of randomness given by the biased WeightDie, it is biased to the implanted motif) although there is a big chance to pick it.
    # This is intentional! It avoids finding a local motif that only looks like being the implanted one
    n = len(Text)
    probabilities = {}
    for index in range(n - k + 1):
        kmer = Text[index: index + k]
        probabilities[kmer] = Pr(kmer, profile)
    return WeightedDie(Normalize(probabilities))

# print(ProfileGeneratedString('AAACCCAAACCC', {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}, 2))


def GibbsSampler(k, t, N, Dna):
    # first, randomly select a k-mer motif from each string in Dna
    BestMotifs = RandomMotifs(Dna, k, t)
    # start loop
    for _ in range(N):
        # then pick one string to delete
        i = randint(1, t) - 1
        deleted_string = Dna[i]
        # then create a profile on the remaining ones
        profile_of_rest = ProfileWithPseudocounts(BestMotifs[:i] + BestMotifs[i + 1:])
        # then extract a motif back from the deleted string plus the profile by using WeightDie randomness
        new_motif = ProfileRandomlyGeneratedString(deleted_string, profile_of_rest, k)
        # recreate the motifs list with the new motif in (the same) place
        renewed_motifs = BestMotifs[:i] + [new_motif] + BestMotifs[i + 1:]
        if Score(BestMotifs) > Score(renewed_motifs):
            # calculate the score and if better than the old ones store the new motifs as best
            BestMotifs = renewed_motifs
    # repeat N times
    # end loop
    return BestMotifs





if __name__ == "__main__":
    import subprocess
    from file_io import outputter, inputter
    with open('../../Downloads/dataset_163_4 (2).txt') as input_file:
        args = [inputter.inputter(word) for line in input_file for word in line.split()]
        # args = [inputter.inputter(line) for line in input_file] # by line

        # index_for_profile = 2
        # profile = {key: list(map(float, args[num + index_for_profile].split(' '))) for (num,key) in enumerate('ACGT') }

    # produce output here
    # output = ProfileMostProbableKmer(args[0].split('\n')[0], args[1], profile)
    i = 0
    BestMotifs = GibbsSampler(*args)
    while i < 20:
        motifs = GibbsSampler(*args)
        if Score(BestMotifs) > Score(motifs):
            BestMotifs = motifs
        i += 1
    output = BestMotifs

    print(outputter.outputter(output))



    with open('output.txt', "w") as output_file:
        output_file.write(outputter.outputter(output))
        # display in default GUI
        subprocess.run(['open', 'output.txt'])




