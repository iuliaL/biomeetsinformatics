# converter from base 10 to some other base
def base_convert(number, base):
    result = []
    while number > 0:
        remainder = number % base
        result.insert(0, str(remainder))  # prepend
        number = number // base  # exclude the remainder
    result_as_string = "".join(result)
    return result_as_string


# print(base_convert(5437, 4))  # => 1110331

# since I have 4 possible nucleotides in DNA
# we will first order all 4k k-mers lexicographically (i.e., according to how they would appear in the dictionary)
# and then convert them into the 4k different integers between 0 and 4k − 1

# A -> 0
# C -> 1
# G -> 2
# T -> 3

# 1110331 -> CCCATTC


def PatternToNumber_(pattern):
    # where number means the index of the pattern in the frequency array
    nucleotides = dict(A=0, C=1, G=2, T=3)
    k = len(pattern)
    i = k - 1  # start from the last nucleotide
    result = 0

    while i >= 0:
        nucleotide = pattern[(k - 1) - i]
        nucleotide_val = nucleotides[nucleotide]
        result += nucleotide_val * (4 ** i)
        i -= 1
    return result

# print(PatternToNumber("ATGCAA")) # -> 912

# How it works?

# ex: PatternToNumber("ATGCAA")
# in base 4 ATGCAA -> 032100

# [0*(4^5)] + [3*(4^4)] + [2*(4^3)] + [1*(4^2)] + [0*(4^1)] + [0*(4^0)] = 0 + 768 + 128 + 16 + 0 + 0 = 912

# or easier using int(string, base)

def PatternToNumber(pattern):
    nucleotides = dict(A=0, C=1, G=2, T=3)
    # write the pattern as in base 4
    as_base_4_string = "".join(
        [str(nucleotides[x]) for x in pattern]
    )
    return int(as_base_4_string, 4)  # parses a number from string-number into the provided base returns base 10

def RecursivePatternToNumber(pattern):
    if len(pattern) == 0:
        return 0
    else:
        nucleotides = dict(A=0, C=1, G=2, T=3)
        return 4 * PatternToNumber(pattern[:-1]) + nucleotides[pattern[-1]]
    
# print(RecursivePatternToNumber("CGACCACCCTATGGCATTC")) # -> 912


# Transform a frequency array index in a pattern
# k is the pattern length

def NumberToPattern(number, k):
    as_base_4_string = base_convert(number, 4)
    length = len(as_base_4_string)
    # take the kth from the right or pad left with zeros
    correct_length_string = as_base_4_string[length - k:] if length > k else as_base_4_string.zfill(k)
    # replace with nucleotides
    nucleotides = ['A', 'C', 'G', 'T']
    pattern = ''
    for string_number in correct_length_string:
        pattern += nucleotides[int(string_number)]
    return pattern


def RecursiveNumberToPattern(number, k):
    nucleotides = ['A', 'C', 'G', 'T']
    quotient = number // 4
    remainder = number % 4
    if k > 1:
        return RecursiveNumberToPattern(quotient, k - 1) + nucleotides[remainder]
    else:
        return nucleotides[remainder]

# print(NumberToPattern(6599, 11))  # -> AAAACGCTACT
# print(RecursiveNumberToPattern(5437, 8))  # -> ACCCATTC

# 4 nucleotides * k positions means 4**k possible kmers out of the 4 nucleotides
def ComputingFrequencies(Text, k):
    # initiate the frequency array with zeros
    frequency_array = [0] * 4**k
    # go through the Text sliding kmer windows till the end and store the occurences of each kmer in the frequency array
    for index in range(len(Text) - k + 1):
        kmer = Text[index: index + k]
        kmer_index = PatternToNumber(kmer)  # => 644
        frequency_array[kmer_index] += 1
    return frequency_array


text = "ATTATTCTTAGGATGAGGCCCGTAGTGCAACGTACATCGCGGGTCAGAGCGATTCGACGTACTTCTGGTCGCGGGGGTGACGCTCCCTGCTTTAACCCGAAAGGAGATCGCGTTTCCTAGCTAAGTCATACAAGTCGGGGCCCTACTCCTCCAGACCCCTAAATCGGACTTGGTCGTCAGTAACCTTAATGCGCTCTTGAACCACTCGCACCTTCCGCATCTGCGGAAGTCCCAGTTCTCGTTGCTTAAGTGGAGGCACACTGCGTGGCCACCATGAAAAGGGTATCTGTTGGACTTTTGGGTCAATTCTATCTGCCCTCGGCACAAAAAAGGAAGTACCCCACTACACCATCGTCTGGTAGGGAACACTGATTTACTTACATAGACCTGCGTCACACTCAACATCTGCCTAAGGAAGGAATTTTGTACGAGACGGTAATTATGATGTTGCATGGCCAGCGGGGCGGGATCTTTGAACATTGTGGCAGGGCAAGGTGCCTCCTATGAGGATCTGCCATCCTTGCTAGCTAGCCGGTATTACGCGCCCGATTAATTTGTCTAAACACAGAATCTTATAAATGAAGACACCGCACTGAGACGGGGGG"
frequency = ComputingFrequencies(text, 6)
output = " ".join(str(x) for x in frequency)

if __name__ == "__main__":
    import subprocess

    file = open('output.txt', "w")
    file.write(output)
    file.close()
    # display in default GUI
    #subprocess.run(['open', 'output.txt'])
