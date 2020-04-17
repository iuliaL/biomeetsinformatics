# converter from base 10 to some other base


def base_convert(number, base):
    result = []
    while number > 0:
        remainder = number % base
        result.insert(0, str(remainder))  # prepend
        number = number // base  # exclude the remainder
    result_as_string = "".join(result)
    return result_as_string


print(base_convert(5437, 4))  # => 1110331

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

# or easier using int(string, base)


def PatternToNumber(pattern):
    nucleotides = dict(A=0, C=1, G=2, T=3)
    # write the pattern as in base 4
    as_base_4_string = "".join(
        [str(nucleotides[x]) for x in pattern]
    )
    return int(as_base_4_string, 4)  # parses a number from string-number into the given base


# print(PatternToNumber("ATGCAA")) # -> 912

# ex: PatternToNumber("ATGCAA")

# in base 4 ATGCAA -> 032100

# [0*(4^5)] + [3*(4^4)] + [2*(4^3)] + [1*(4^2)] + [0*(4^1)] + [0*(4^0)] = 0 + 768 + 128 + 16 + 0 + 0 = 912

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


print(NumberToPattern(5437, 7))  # -> CCCATTC
print(NumberToPattern(5437, 8))  # -> ACCCATTC
