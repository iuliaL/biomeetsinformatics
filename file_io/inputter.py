# parse string arguments as float, int or just string
def inputter(arg):
    if '.' in arg:
        try:
            return float(arg)
        except ValueError:
            return arg
    else:
        try:
            return int(arg)
        except ValueError:
            return arg