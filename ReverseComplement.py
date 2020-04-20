# Reverse Complement Problem: Find the reverse complement of a DNA string.
#      Input: A DNA string Pattern.
#      Output: The reverse complement of Pattern.


def ReverseComplement(Text):
    pairs = dict(A="T", C="G", G="C", T="A")
    complement = [
        pairs[x] for x in Text
    ]
    return "".join(complement[::-1])


output = ReverseComplement("ATGATCAAG")  # => CTTGATCAT

if __name__ == "__main__":
    import subprocess

    file = open('output.txt', "w")
    file.write(output)
    file.close()
    # display in default GUI
    subprocess.run(['open', 'output.txt'])
