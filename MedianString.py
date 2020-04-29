from HammingDistance import HammingDistance
from Combinations import Combinations

# this function might not be right
def MedianString(k, *Dna):
    median = None
    for comb in Combinations(k):
        for string in Dna:
            distance = float('inf')
            for i in range(len(string) - k + 1):
                curr_kmer = string[i: i + k]
                curr_distance = HammingDistance(curr_kmer, comb)
                if distance > curr_distance:
                    median = comb
                    distance = curr_distance
    return median


if __name__ == "__main__":
    import subprocess
    from outputter import outputter
    from inputter import inputter
    with open('../../Downloads/dataset_158_9 (4).txt') as input_file:
        args = [inputter(word) for line in input_file for word in line.split()]

    # produce output here
    output = MedianString(*args)

    with open('output.txt', "w") as output_file:
        output_file.write(outputter(output))

    # display in default GUI
    subprocess.run(['open', 'output.txt'])

