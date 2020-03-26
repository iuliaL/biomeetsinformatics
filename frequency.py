# k is the length of a k-mer(Pattern)

def FrequencyMap(Text, k):  
    frequency = {}
    for index in range(len(Text) - k + 1):
        Pattern = Text[index: index + k]
        # print("k-mer is " + Pattern)
        if Pattern not in frequency:
            frequency[Pattern] = 1
        else:
            frequency[Pattern] += 1
    return frequency

# ( ^ the well known occurences hashmap problem )

def FrequentPatterns(Text, k):
    most_frequent_patterns = []
    frequency = FrequencyMap(Text, k)
    at_least_3_times = getFrequent(frequency, 3)
    print(at_least_3_times)
    maximum_frequency = max(frequency.values())
    print("How frequent is max? " + str(maximum_frequency))

    for pattern in frequency:
        if frequency[pattern] == maximum_frequency:
            most_frequent_patterns.append(pattern)
    return most_frequent_patterns

def getFrequent(frequencies, threshold):
    frequent = {}
    for key in frequencies:
        if frequencies[key] >= threshold:
            frequent[key] = frequencies[key]
    return frequent

""" 
Vibrio Cholerae oriC

ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCA
CACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCC
TTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGG
CCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATAT
TGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTT
TATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGC
CTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTAC
CTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCA
TGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC

"""

vibrio_cholerae_oriC_file = open(r"vibrio_cholerae_oriC.txt")
vibrio_cholerae_oriC = vibrio_cholerae_oriC_file.read()
vibrio_cholerae_oriC_file.close()


# print(FrequentPatterns(vibrio_cholerae_oriC, 9)) # => ['ATGATCAAG', 'CTCTTGATC', 'TCTTGATCA', 'CTTGATCAT']
""" Note the first 2 are complementary strands """


t_petrophila_oriC_file = open(r"t_petrophila_oriC.txt")
t_petrophila_oriC = t_petrophila_oriC_file.read()
t_petrophila_oriC_file.close()
print(FrequentPatterns(t_petrophila_oriC, 9))