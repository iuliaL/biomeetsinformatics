def outputter(output):
    if  isinstance(output, list):
        return "\n".join([
            str(x) for x in output
        ])
    else:
        return str(output)
