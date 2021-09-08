# created at sep 8th to extract variant id from uniprot txt API, store and insert in expasy API to
# get the mutation position fromn there
import urllib.request

url = 'https://www.uniprot.org/uniprot/Q14204.txt'
file_test = urllib.request.urlopen(url)
for l in file_test:  # TODO: include a condition that first checks if the DISEASE lines exist, decide what to store from that part, check how many of your uniprot IDs remain based on this condition
	line = l.decode("utf-8")
	if line.startswith('CC') and '-!- DISEASE:' in line:
		disease_name = line.split('CC   -!- DISEASE:')[1]
		if '[' in line:
			disease_name = disease_name.split('[')[0]
			break
		else:
			line = line.next
			disease_name = disease_name.append()
		print(disease_name)
		# if line.startswith('FT') and 'VARIANT' in line:
		# 	print('There is a variant here!')
		# 	position = line.split('VARIANT')[1:][0].replace(' ', '')
		# 	print(position)

	# print(decoded_line)