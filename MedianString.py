from HammingDistance import HammingDistance
from Combinations import Combinations

def MedianString(k, *Dna):
    median = None
    distance = float('inf') # infininty
    for comb in Combinations(k):
        curr_distance = DistanceBetweenPatternAndAllDnaStrings(comb, *Dna)
        if distance > curr_distance:
            distance = curr_distance
            median = comb
    return median

def DistanceBetweenPatternAndAllDnaStrings(Pattern, *Dna):
    k = len(Pattern)
    total_distance = 0
    for string in Dna:
        distance = float('inf')
        for i in range(len(string) - k + 1):
            curr_kmer = string[i: i + k]
            curr_distance = HammingDistance(curr_kmer, Pattern)
            if distance > curr_distance:
                distance = curr_distance
        total_distance += distance
    return total_distance

# if __name__ == "__main__":
#     import subprocess
#     from file_io import outputter, inputter
#     with open('../../Downloads/dataset_158_9 (5).txt') as input_file:
#         args = [inputter.inputter(word) for line in input_file for word in line.split()]

#     # produce output here
#     output = MedianString(*args)

#     with open('output.txt', "w") as output_file:
#         output_file.write(outputter.outputter(output))
#         subprocess.run(['open', 'output.txt']) # display in default GUI




