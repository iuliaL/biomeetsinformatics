import time
from FrequentWords import FrequencyMap, getFrequent

# k = length of kmer (pattern)
# t = times occuring
# L =  window length (what defines a clump)


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
        last_first_kmer = genome[index - 1 : index - 1 + k]

        new_last_kmer = genome[index + L - k : index + L]
        # decrease the last kmer freq by 1
        frequencies_in_window[last_first_kmer] -= 1
        # increase the new kmer freq by 1 if exists else set it to 1
        if new_last_kmer not in frequencies_in_window:
            frequencies_in_window[new_last_kmer] = 1
        else:
            frequencies_in_window[new_last_kmer] += 1
        
        if frequencies_in_window[new_last_kmer] >= t:
            kmers_forming_clumps.add(new_last_kmer)


    end = time.process_time()
    print("Elapsed time:", end - start)
    return len(list(kmers_forming_clumps))


# genome = "CCATAAGCACATCCCCTCAGTCAGGTATGGCACTTAACCGGTGTAGGCTTGGATGGCCCAGGTGCTTGTGCTATCCGTAGCCGTGGTCAATCATCTGCCATAAAACATAATGTCTGATTAGTACGGACTTTTGTGGTGATGACTATCATCGTGAGACCGTTCTATCCTCCAGTCGTACCCGGTAATCGCTCCTTTACGCGATTGTACTTCGGCTCCCGGCAATAACAAGTAGTAGGAAGAAACCCTGTAATAGTGAAACAGTGTCCTCAAGCAAACTACATCCAACGGTCAACTGTATCTACGAGATCTCGAAGGAGGTAGGTGGACCGCCGTGGTGCTGTCCGTTTGCACAGGGGATGACAGTCAACCGCGAATGGAGCAGGGTGCTCGATTGTTACACTCGGTGTACTTCTGAAAAGACCGCAAACTGTTCAATATCTGGGCATGATGCATGATGCATGATCCAAACAGGGACCCACATATGCTGATCGTGTCAGTCCATATCTCACAACGTCATTAGTTGCCATGATAATGTTCCGCCTTGTAACTCTCTGGTTGCCTCATCCTCTTGGAATGAGTCGTGACGCACACCAATACTGGATTTTTCACTTTTTCACTTTTTTCACTTGTCGCACCGCTAAGCCTAGGAGTTCGTGAGCTGGCTGTGTCACCTGCGCCAGGTGATTGCTTCCGCAAATATCGGAACAATGAGGTGTGACGCTGTAGGGGACCCCTTTTAGTGGATGGAGATCTATTAGTAGGACAGCCCTCCGGAAATCACTATACGATGAGCCCGATGAAGTGTTGATAAGCCCCCGCGGAATTCGGCCTATTTATCCGGGTATCCTCAATCACGTGTGAGACTGATTTCGTCATCTATGAAGGGAGGCCTGCAGTTGAAATTCTTCTTGAGGGTCATCCCCCAACACTCTCCACCAACACTCCAACACTTCGTCGTCCGAAAGAGACAGAGGACCCTACAGTTAATGCCCAAGAATCTTAACATTGTTCTGTCATTGTTCCCCATTGTTCCAATACTTAAAGGCGCCATAGCATTACCAGTTCAGACGGGGTGTGTATCTCCACTCGAATCCTGCCACATACAGATACATACAGATTGAAGTAGGATGGATTTAGGATATGTGTTGTTGTCGCGGGGCGGCGCGGGGTTGAGTGTACTGCGTGTCTCCCCAAGAGGCAGTACAACAGTCAATTGCAAGAGGCCAGAAGTTGCGGACCATCATGTCTAATGCATAACGCTGCGCGCTATCTTTCTCTAGAAGAGAGTTATCTGTATTGTTGACTCACCACCTATGGGTTCCCCCCGGTTCCCCGGGTTCCCCGGAGGGCTGGCGCACTTCTATAGGACACTTGGGGATTTGACCTGACAGATGTAGCATTTCGACATCATCTCACATCTTTCTCTTTCTTCAACGTTCAGTCAGCTCAAGCTCAAGCTCATTACGAACTGTAGATCGTAAAGTCATGTTTATAAGCAGGAAGAAATTTACACCACTAAATTGCTCCGCGCGGGTGCGCGGGTGCGCGGGTGGGTGCATAGAGTCTTCCTTACCCTAACTAGCTTAGCTGCAGTTAACCACTTCAAGGGAACGAGACTCGCATCCCTTTGCGCTGTGTGTGAGGTACGGGGGGGGGGGGGGAGGGGGGGAGAGGAGTTGTCCGACTGGGGCAGTGTCCAGCAGACCGTCAAAGCAAGATATTCAGAAGATCAGAAGAAAGTCAGAAGAGCAACTTCTTGTATTTAAGACAGGCTTAAGTACGGGTGCTGATCTAAATGTCTGAGCCTAACCAGGATGGTCGCGATTTGACTAACTAAGCAATATCGTCACCATTCAACTACGCTATAGCTTGCGAGCTTGCGTTGCGGCCGTTTAGCGATCCACGGTAGCACGGTAGCACGGTAGCACGGTAGCGTAACGGGGTAACGGGGTAACGGGGTAACGGGGTAACGGGGTAACGGGGTAACGGTAACGGG"
# problem = ClumpFinding(genome, 8, 25, 4)
# output = " ".join(problem)
# print(output)

# if __name__ == "__main__":
#     import subprocess

#     file = open('output.txt', "w")
#     file.write(output)
#     file.close()
#     # display in default GUI
#     subprocess.run(['open', 'output.txt'])

# Problem:
# How many different 9-mers form (500,3)-clumps in the E. coli genome? (In other words, do not count a 9-mer more than once.)

with open('genome/E_coli_genome.txt') as file:
    e_coli = file.read()
    print("Answer:", ClumpFinding(genome=e_coli, k=9, L=500, t=3))

# Elapsed time: 10.338602
# Answer: 1904