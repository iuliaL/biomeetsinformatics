import time
from collections import defaultdict

# k is the length of a k-mer(Pattern)


def FrequencyMap(Text, k):
    frequency = defaultdict(lambda: 0)
    for index in range(len(Text) - k + 1):  # exp: "ATCCGA" in 6 length text i have 4 3-length kmers 6 - 3 + 1 = 4; range 0 -> 4 means 4 times
        Pattern = Text[index: index + k]
        frequency[Pattern] += 1
    return frequency

# ( ^ the well known occurences hashmap problem )


# filter only the results that have maximum frequency
def FrequentWords(Text, k):
    # start = time.clock()
    most_frequent_patterns = []
    frequency = FrequencyMap(Text, k)
    maximum_frequency = max(frequency.values())
    # print("How frequent is max? " + str(maximum_frequency))

    for pattern in frequency:
        if frequency[pattern] == maximum_frequency:
            most_frequent_patterns.append(pattern)
    # end = time.clock()
    # print("Elapsed time:", end - start)
    return most_frequent_patterns

# print(FrequentWords("CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA", 3))


def getFrequent(frequencies, threshold):
    # filter only the results appearing more often than threshold
    frequent = []
    for key in frequencies:
        if frequencies[key] >= threshold:
            frequent.append(key)
    return frequent

