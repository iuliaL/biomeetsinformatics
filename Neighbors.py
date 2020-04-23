import itertools


def ImmediateNeighbors(pattern):
    """
    Gets all strings within hamming distance 1 of motif and returns it as a
    list.
    """
    immediate_neighbors = set()
    nucleotides = ["A", "C", "G", "T"]
    for index in range(len(pattern)):
        nucleotide = pattern[index]
        for n in nucleotides:
            if not nucleotide == n:
                neighbor = pattern[:index] + n + pattern[index + 1:]
                immediate_neighbors.add(neighbor)
    return immediate_neighbors


# print(ImmediateNeighbors("TGCA"))

def Neighbors(motif, d):  # where d is hamming distance
    workingSet = {motif}
    for _ in range(d):
        workingSet = set(itertools.chain.from_iterable(map(ImmediateNeighbors, workingSet)))
    return list(workingSet)

