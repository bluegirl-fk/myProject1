# created at sep 8th to extract variant id from uniprot txt API, store and insert in expasy API to
# get the mutation position fromn there
import urllib.request

url = 'https://www.uniprot.org/uniprot/Q14204.txt'
file_test = urllib.request.urlopen(url)
for l in file_test:
	line = l.decode("utf-8")
	if line.startswith('FT') and 'VARIANT' in line:
		print('There is a variant here!')
	# print(decoded_line)