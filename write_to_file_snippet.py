output = ""

if __name__ == "__main__":
    import subprocess

    file = open('output.txt', "w")
    file.write(output)
    file.close()
    # display in default GUI
    subprocess.run(['open', 'output.txt'])