# Entropy is a measure of the uncertainty of a probability distribution (p1, . . . , pN)
 
# For example, the entropy of the probability distribution (0.2, 0.6, 0.0, 0.2)

# -(0.2 \log_{2}{0.2} + 0.6\log_{2}{0.6} + 0.0\log_{2}{0.0} + 0.2\log_{2}{0.2}) \approx 1.371
from math import log2

def Entropy(*probabilities):
    entropy = 0
    for p in probabilities:
        if p > 0:
            entropy += p * log2(p)
    return entropy * -1

print(Entropy(0.2, 0.6, 0.0, 0.2)) # => 1.3709505944546687

a=[0.2,0.2,0.9,0.1,0.1,0.1,0.3]
c=[0.1,0.6,0.4,0.1,0.2,0.4,0.6]
g=[1,1,0.9,0.9,0.1]
t=[0.7,0.2,0.1,0.1,0.5,0.8,0.7,0.3,0.4]

 # collective entropy
print(Entropy(*a,*c,*g,*t))
        
