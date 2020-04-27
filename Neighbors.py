from HammingDistance import HammingDistance

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

def Neighbors(pattern, d):
    neighbors =  set()
    neighbors.add(pattern) # initialize neighbors with pattern
    for _ in range(d):
        for n in neighbors:
            neighbors = neighbors | ImmediateNeighbors(n) # (set union)
    return list(neighbors)

# HORRIBLE
def RecursiveNeighbors(pattern, d):
    k = len(pattern)
    nucleotides = ["A", "C", "G", "T"]
    if d == 0:
        return { pattern }
    if k == 1 :
        return set(nucleotides)

    neighbors = set()
    suffix = pattern[1:]
    first_symbol = pattern[0]
    suffix_neighbors = RecursiveNeighbors(suffix, d)
    for s_neighbor in suffix_neighbors:
        distance = HammingDistance(suffix, s_neighbor)
        if distance < d:
            for n in nucleotides:
                neighbor = n + s_neighbor
                neighbors.add(neighbor)
        else:
            neighbor = first_symbol + s_neighbor
            neighbors.add(neighbor)
    return list(neighbors)


