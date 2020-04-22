

if __name__ == "__main__":
    import subprocess
    from outputter import outputter
    from inputter import inputter
    with open('../../Downloads/dataset_9_3.txt') as input_file:
        args = [inputter(word) for line in input_file for word in line.split()]

    # produce output here
    output = "_______Missing output!______"

    with open('output.txt', "w") as output_file:
        output_file.write(outputter(output))

    # display in default GUI
    subprocess.run(['open', 'output.txt'])


