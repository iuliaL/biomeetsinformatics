#  Hamming Distance Problem:
#  The total number of mismatches between strings p and q is called the Hamming distance between these strings.
#  Compute the Hamming distance between two equal length strings.
#  Input: Two strings of EQUAL length.
#  Output: The Hamming distance between these strings.



def HammingDistance(p, q):
    mismatches = 0
    for i in range(len(p)):
        if not (p[i] == q[i]):
            mismatches +=  1
    return mismatches 

