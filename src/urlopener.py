# created at sep 8th to extract variant id from uniprot txt API, store and insert in expasy API to
# get the mutation position fromn there
import urllib.request

url = 'https://www.uniprot.org/uniprot/Q14204.txt'
file_test = urllib.request.urlopen(url)
for l in file_test:  # TODO: include a condition that first checks if the DISEASE lines exist, decide what to store from that part
	line = l.decode("utf-8")
	if line.startswith('FT') and 'VARIANT' in line:
		print('There is a variant here!')
		position = line.split('VARIANT')[1:][0].replace(' ', '')
		print(position)

	# print(decoded_line)