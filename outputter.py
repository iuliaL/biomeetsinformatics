def outputter(output):
    if  isinstance(output, list):
        return " ".join([
            str(x) for x in output
        ])
    else:
        return str(output)
