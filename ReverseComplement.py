# Reverse Complement Problem: Find the reverse complement of a DNA string.
#      Input: A DNA string Pattern.
#      Output: The reverse complement of Pattern.


def Complement(text):
    rule = {
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C"
    }
    new_strand = ""
    for nucleotide in text:
        new_strand += rule[nucleotide]
    return new_strand

def Reverse(text):
    # slice with negative step goes backwards
    return text[::-1]

def ReverseComplement(text):
    return Reverse(Complement(text))


print("ATGATCAAG has the complement strand " + ReverseComplement("ATGATCAAG")) # => CTTGATCAT

