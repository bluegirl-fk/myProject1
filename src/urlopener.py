# created at sep 8th to extract variant id from uniprot txt API, store and insert in expasy API to
# get the mutation position fromn there
import urllib.request
import re


url = 'https://www.uniprot.org/uniprot/Q14204.txt'
my_par = []
file_test = urllib.request.urlopen(url)
for l in file_test:
	line = l.decode("utf-8")
	if re.search(r'^CC\s+-!-', line):
		if re.search(r'^CC\s+-!- DISEASE:', line):
			# Crerate new paragraph
			my_par.append(line)
		# Otherwise, break only if all DISEASE line have been terminated
		elif len(my_par) > 0:
			break
	elif len(my_par) > 0:
		my_par[-1] += line

print(my_par)

