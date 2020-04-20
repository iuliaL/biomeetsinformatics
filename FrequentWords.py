import time

# k is the length of a k-mer(Pattern)

def FrequencyMap(Text, k):
    frequency = {}
    for index in range(len(Text) - k + 1): # exp: "ATCCGA" in 6 length text i have 4 3-length kmers 6 - 3 + 1 = 4; range 0 -> 4 means 4 times
        Pattern = Text[index: index + k]
        # print("k-mer is " + Pattern)
        if Pattern not in frequency:
            frequency[Pattern] = 1
        else:
            frequency[Pattern] += 1
    return frequency

# ( ^ the well known occurences hashmap problem )


# filter only the results that have maximum frequency
def FrequentWords(Text, k):
    # start = time.clock()
    most_frequent_patterns = []
    frequency = FrequencyMap(Text, k)
    # at_least_3_times = getFrequent(frequency, 3)
    # print(at_least_3_times)
    maximum_frequency = max(frequency.values())
    # print("How frequent is max? " + str(maximum_frequency))

    for pattern in frequency:
        if frequency[pattern] == maximum_frequency:
            most_frequent_patterns.append(pattern)
    # end = time.clock()
    # print("Elapsed time:", end - start)
    return most_frequent_patterns


def getFrequent(frequencies, threshold):
    # filter only the results appearing more often than threshold
    frequent = []
    for key in frequencies:
        if frequencies[key] >= threshold:
            frequent.append(key)
    return frequent


# if __name__ == "__main__":

#     vibrio_cholerae_oriC_file = open(r"oriC/vibrio_cholerae_oriC.txt")
#     vibrio_cholerae_oriC = vibrio_cholerae_oriC_file.read()
#     vibrio_cholerae_oriC_file.close()


#     # print(FrequentWords(vibrio_cholerae_oriC, 9)) # => ['ATGATCAAG', 'CTCTTGATC', 'TCTTGATCA', 'CTTGATCAT']
#     """ Note the first 2 are complementary strands """


#     t_petrophila_oriC_file = open(r"oriC/t_petrophila_oriC.txt")
#     t_petrophila_oriC = t_petrophila_oriC_file.read()
#     t_petrophila_oriC_file.close()
#     print(FrequentWords(t_petrophila_oriC, 9))
