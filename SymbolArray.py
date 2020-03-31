# Find the occurences of given nucleotide "symbol" (A, T, C, G) in a window (genome substring) of length equal to half of the genome
# the point is to compare the nucleotide occurences in the 2 halfs of the genome in order to find the oriC.
# The point where the C starts increasing after having decreased (see plot) might indicate the OriC
# (remember! the Genome of Bacteria is CIRCULAR so the index = 0 does not mean the genome starts at that position 0! )
# we have to extend the genome be a half to mimic the "circularity"
# return an array (like an indexed dict)

from PatternCount import PatternCount
import matplotlib.pyplot as plt

def SymbolArray(Genome, symbol):
    n = len(Genome)
    # window_length = n // 2 # => integer division (get rid of remainder)
    # window = Genome[0 : n // 2]
    ExtendedGenome = Genome + Genome[0:n//2]
    symbolArray = {}
    symbolArray[0] = PatternCount(Genome[0:n//2], symbol)
    for index in range(1, len(Genome)):
        found_symbol_in_previous_1st_position = - 1 if symbol == ExtendedGenome[index - 1] else 0 
        found_symbol_in_curr_last_position = 1 if symbol == ExtendedGenome[index + n // 2 - 1] else 0
        occurences = symbolArray[index - 1] + found_symbol_in_previous_1st_position + found_symbol_in_curr_last_position
        symbolArray[index] = occurences
    return symbolArray

with open('E_coli.txt') as file:
    e_coli = file.read()

array = SymbolArray(e_coli, "C")
max_occurences = max(array.values())
print("Max C count:", max_occurences) # => 606875

# plt.plot(*zip(*sorted(array.items()))) # don't know why sorted since this is actually a list
plt.plot(array.values())
plt.xlabel('Genome position')
plt.ylabel('C count in half-genome legth window')
plt.show()

if __name__ == "__main__":
    import sys

    with open(sys.argv[1],'r') as f:
        genome = f.read()
    symbol =  sys.argv[2]
    SymbolArray(genome, symbol)