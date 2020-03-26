# import re

def PatternCount(Text, Pattern):
	count = 0
	i = 0
	while i < len(Text) - len(Pattern) + 1:
		curr =  Text[i: i + len(Pattern)]
		# print('Curr   ' + str(i) + " " + curr)
		if curr  == Pattern:
			count += 1
		i += 1
	return count
	



pattern = "ATA"
text = "CGTAATATCCATAG"

count = PatternCount(text, pattern)
print("Pattern found: {}".format(count)) # => 2

# def PatternCountWithRegex (Text, Pattern):
# 	return len(re.findall(Pattern, Text))

# count = PatternCountWithRegex(text, pattern)

# print("Pattern found: {}".format(count))