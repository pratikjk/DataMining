import csv
data_csv_file="GeneAssocaitionData.csv"
gene_prefix="G"

def one_item_frequentSet(data_set):
	item_count=dict()
	for row in data_set:
		for item in row:
			if item not in item_count:
				item_count[item] = 1
			else:
				item_count[item] += 1
	return item_count

def apriori_algorithm(data_set):
	# get frequent itemset of 1 items
	C=list()
	L=list()
	C.append(one_item_frequentSet(data_set))
	L.append(prune(C[0]))
	k=1
	
	while not C[k-1]:
		C[k].append(join(L[k-1],L[k-1]))
		L[k].append(prune(C[k]))
		k=k+1
	
	return C



# @task: add support for reading filename from commandline input
try:
	data_file=open(data_csv_file,"rU")
	try:
		reader = csv.reader(data_file)
		modified_data=[]
		for row in reader:
			colnum=1
			modified_row=[]
			for col in row:
				modified_row.append(gene_prefix+str(colnum)+col)
				colnum=colnum+1
			modified_data.append(modified_row)
		#print (modified_data)
		#line=[['pratik','kathalkar','computer'],['computer','science','engineering'],['pratik','coep','engineering','and','science','college']]
		print(len(one_item_frequentSet(modified_data)))
	finally:
		data_file.close()
except IOError:
	print('An error occured trying to read the file, file may not be present')



