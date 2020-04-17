# OFF TOPIC

# Summing Integers (to n) Problem
def RecursiveGauss(n):
    if n == 0:
        return 0
    else:
        return RecursiveGauss(n-1) + n

print(RecursiveGauss(5))

def GeniusGauss(n):
    return n * (n + 1) // 2 # (I can use integer division since n is integer and either n or n+1 is even)

print(GeniusGauss(5))





