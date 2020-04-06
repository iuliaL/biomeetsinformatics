# Input:
    # AACGTA
    # CCCGTT
    # CACCTT
    # GGATTA
    # TTCCGG

from HammingDistance import HammingDistance

def CountMatrix(Motifs):
    count = {}
    k = len(Motifs[0])
    for key in "ACGT":
        count[key] = [ 0 ] * k
    for motif in Motifs:
        for index in range(k):
            count[motif[index]][index] += 1
    return count


print(CountMatrix(['AACGTA','CCCGTT', 'CACCTT', 'GGATTA','TTCCGG']))


def Profile(Motifs):
    count = CountMatrix(Motifs)
    for nucleotide_key in count:
        count[nucleotide_key] = [x / len(Motifs) for x in count[nucleotide_key]]
    return count

# print(Profile(['AACGTA','CCCGTT', 'CACCTT', 'GGATTA','TTCCGG'] ))

def Consensus(Motifs):
    consensus = ''
    count = CountMatrix(Motifs)
    k = len(Motifs[0])
    for index in range(k):
        max_so_far = 0
        most_freq_nucleotide = None
        for nucleotide, val in count.items():
            if val[index] > max_so_far:
                max_so_far = val[index]
                most_freq_nucleotide = nucleotide
        consensus += most_freq_nucleotide
    return consensus

print(Consensus(['AACGTA','CCCGTT', 'CACCTT', 'GGATTA','TTCCGG'])) # => CACCTA

def Score_(Motifs):
    consensus = Consensus(Motifs)
    count = CountMatrix(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    score = 0
    for index in range(k):
        subtract = 0
        nucleotide = consensus[index]
        for key,val in count.items():
            if key == nucleotide:
                subtract += val[index]
        score += t - subtract
    return score

def Score(Motifs):
    score = 0
    consensus = Consensus(Motifs)
    for motif in Motifs:
        score += HammingDistance(motif, consensus)
    return score
        

print(Score(['AACGTA','CCCGTT', 'CACCTT', 'GGATTA','TTCCGG']))
    
    
print(Score(["GTACAACTGT",
"CAACTATGAA",
'TCCTACAGGA',
'AAGCAAGGGT',
'GCGTACGACC',
'TCGTCAGCGT',
'AACAAGGTCA',
'CTCAGGCGTC',
"GGATCCAGGT",
"GGCAAGTACC"]))
            
            
        
        

            