import time
from FrequentWords import FrequencyMap, getFrequent
from collections import defaultdict

# k = length of kmer (pattern)
# t = times occuring
# L =  window length (what defines a clump)

# Problem
# Count how many different kmers appear in a window length L at least t times (clump a region where kmers repeat often)

# def ClumpFinding_(genome, k, L, t):
#     start = time.process_time()
#     kmers_forming_clumps = []
#     for index in range(len(genome) - L + 1):
#         window = genome[index: index + L]
#         kmer_frequency_in_window = FrequencyMap(window, k)
#         at_least_t_times = getFrequent(kmer_frequency_in_window, t)
#         kmers_forming_clumps += at_least_t_times

#     de_dupe_kmers = list(set(kmers_forming_clumps))
#     end = time.process_time()
#     print("Elapsed time:", end - start)
#     return de_dupe_kmers

def ClumpFinding(genome, k, L, t):
    start = time.process_time()

    first_window = genome[0:L]
    frequencies_in_window = FrequencyMap(first_window, k) # kmer frequencies in 1st window
    at_least_t_times = getFrequent(frequencies_in_window, t)
    kmers_forming_clumps = set(at_least_t_times)

    for index in range(1, len(genome) - L + 1):
        # slide the window and decrease the last kmer frequency by 1
        last_first_kmer = genome[index - 1 : index - 1 + k]
        frequencies_in_window[last_first_kmer] -= 1

        new_last_kmer = genome[index + L - k : index + L]
        # increase the new kmer freq by 1 if exists else set it to 1
        if new_last_kmer not in frequencies_in_window:
            frequencies_in_window[new_last_kmer] = 1
        else:
            frequencies_in_window[new_last_kmer] += 1
        
        if frequencies_in_window[new_last_kmer] >= t:
            kmers_forming_clumps.add(new_last_kmer)

    end = time.process_time()
    print("Elapsed time:", end - start)

    return len(kmers_forming_clumps)

def FasterClumpFinding(genome, k, L, t):
    start = time.process_time()
    count = 0
    # 1. build dictionary mapping all possible kmers to their indices (3s total).
    indexed_kmers = defaultdict(lambda: [])
    for index in range(len(genome) - k + 1):
        indexed_kmers[genome[index: index + k]].append(index)

    # 2. increment count when the first t indices of a given kmer fit within window length L (3s total).
    for _, index_list in indexed_kmers.items():
        if len(index_list) >= t:
            i = 0
            while i < len(index_list) - t + 1:
                if index_list[i + t - 1] - index_list[i] < L - k + 1:
                    count += 1
                    break
                i+=1
           
    # 3.  return count
    end = time.process_time()
    print("Elapsed time", end - start)
    return count




# Problem:
# How many different 9-mers form (500,3)-clumps in the E. coli genome? (In other words, do not count a 9-mer more than once.)

with open('genome/E_coli_genome.txt') as file:
    e_coli = file.read()
    print("Answer ClumpFinding", ClumpFinding(genome=e_coli, k=9, L=500, t=3))
    print("Answer FasterClumpFinding", FasterClumpFinding(genome=e_coli, k=9, L=500, t=3))


# using ClumpFinding
# Elapsed time: 10.338602
# Answer: 1904

# using FasterClumpFinding
# Elapsed time: 9.726102
# Answer: 1904