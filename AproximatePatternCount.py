# Count how many times Pattern and slightly different Pattern (with d allowed differences) are found in Genome

from HammingDistance import HammingDistance

def ApproximatePatternCount(Pattern, Genome, d):
    k = len(Pattern)
    count = 0
    for index in range(len(Genome) - k + 1):
        curr_kmer = Genome[index:index + k]
        if HammingDistance(Pattern, curr_kmer) <= d:
            count += 1 
    return count 


# print(ApproximatePatternCount('ACAA', 'AACAAGCTGATAAACATTTAAAGAG', 1))  # => 4

# if __name__ == "__main__":
#     import subprocess
#     from outputter import outputter
#     from inputter import inputter
#     with open('../../Downloads/dataset_9_6.txt') as input_file:
#         args = [inputter(word) for line in input_file for word in line.split()]

#     # produce output here
#     output = ApproximatePatternCount(*args)

#     with open('output.txt', "w") as output_file:
#         output_file.write(outputter(output))

#     # display in default GUI
#     subprocess.run(['open', 'output.txt'])