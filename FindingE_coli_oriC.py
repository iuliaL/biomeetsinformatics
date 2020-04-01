
# Identify the most frequent 9-mers (with 1 mismatch) within a window of length 500 starting at position 3923620 of the E. coli genome
# d = 1


# from AproximatePatternMatching import ApproximatePatternMatching
from AproximatePatternCount import ApproximatePatternCount

def MostFreqAproximate9mers(Genome, start_pos=3923620, window_length=500, d=1):
    k = 9
    window = Genome[start_pos: start_pos + window_length + 1]
    results = {}
    for index in range(window_length):
        current9mer = Genome[index:index + k]
        count = ApproximatePatternCount(current9mer, window, d)
        if count > 0 :
            results[current9mer] = count
        
    print(results)

# TODO: not yet done

if __name__ == "__main__":
    import sys
    if len(sys.argv[1:]) != 1:
        print("You must pass a genome string argument")
        exit()
    with open(sys.argv[1],'r') as file:
        genome = file.read()

        MostFreqAproximate9mers(genome)
