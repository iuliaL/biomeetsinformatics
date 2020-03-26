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
    maximum_frequency = max(frequency.values())
    print("How frequent is max? " + str(maximum_frequency))
    for pattern in frequency:
        if frequency[pattern] == maximum_frequency:
            most_frequent_patterns.append(pattern)
    return most_frequent_patterns

""" Vibrio Cholerae oriC """

vibrio_cholerae = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCA"\
"CACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCC"\
"TTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGG"\
"CCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATAT"\
"TGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTT"\
"TATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGC"\
"CTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTAC"\
"CTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCA"\
"TGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"


print(FrequentPatterns(vibrio_cholerae, 9)) # => ['ATGATCAAG', 'CTCTTGATC', 'TCTTGATCA', 'CTTGATCAT']
""" Note the first 2 are complementary strands """