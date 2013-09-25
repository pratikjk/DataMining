import csv
data_csv_file="GeneAssocaitionData.csv"
gene_prefix="G"


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
		print len(modified_data)
	finally:
		data_file.close()
except IOError:
	print('An error occured trying to read the file, file may not be present')
