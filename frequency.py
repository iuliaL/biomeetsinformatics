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
   

def FrequentWords(Text, k):
    most_frequent_patterns = []
    frequency = FrequencyMap(Text, k)
    maximum_frequency = max(frequency.values())
    for pattern in frequency:
        if frequency[pattern] == maximum_frequency:
            most_frequent_patterns.append(pattern)
    return most_frequent_patterns


print(FrequencyMap("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4))
print(FrequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4))