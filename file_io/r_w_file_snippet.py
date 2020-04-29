if __name__ == "__main__":
    import subprocess
    from file_io import outputter, inputter
    with open('../../Downloads/dataset_9_3.txt') as input_file:
        args = [inputter.inputter(word) for line in input_file for word in line.split()]

    # produce output here
    output = "_______Missing output!______"

    with open('output.txt', "w") as output_file:
        output_file.write(outputter.outputter(output))
        subprocess.run(['open', 'output.txt']) # display in default GUI


