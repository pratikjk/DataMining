import csv
import copy
data_csv_file="GeneAssocaitionData.csv"
gene_prefix="G"
min_support=0.50

# get frequent itemset of 1 items
def one_item_frequentSet(data_set,min_support):
	item_count=dict()
	item_set=set()
	for row in data_set:
		for item in row:
			if item not in item_count:
				item_count[item] = 1
			else:
				item_count[item] += 1
	for key in item_count.keys():
		if (float(item_count[key])/len(data_set)) >= min_support:
			item_set.add(frozenset([key]))

	return item_set

def join(frequent_set1, frequent_set2,level):
	candidate_set=set()
	for itemset1 in frequent_set1:
		for itemset2 in frequent_set2:
			if len(itemset1.union(itemset2))==level:
				candidate_set.add(itemset1.union(itemset2))
	return candidate_set

def prune(data_set, candidate_set,min_support):
	pruned_set=set()
	for item_set in candidate_set:
		count=0
		for row in data_set:
			if item_set.issubset(row):
				count+=1
		if (float(count)/len(data_set)>=min_support):
			pruned_set.add(item_set)
	return pruned_set
		
	

def apriori_algorithm(data_set, min_support):
	
	L=one_item_frequentSet(data_set,min_support)
	k=1
	temp_candidate_set=set()
	temp_frequent_set=L
	while temp_frequent_set != set([]):
		temp_candidate_set=join(temp_frequent_set,temp_frequent_set,k)
		temp_frequent_set=prune(data_set,temp_candidate_set,min_support)
		L.update(temp_frequent_set)
		k=k+1	
	
	return L

def generate_all_subsets(elements,index):
	all_subsets=[]
	if(len(elements) == index):
		all_subsets.append([])
	else:
		all_subsets=generate_all_subsets(elements,index+1)
		item=elements[index]
		temp_subsets=copy.deepcopy(all_subsets)
		for tuples in temp_subsets:
			tuples.append(item)
		all_subsets.extend(temp_subsets)
	return all_subsets

def association_rules():
	pass
	# convert the frequent itemsets to list
	

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
		db=[['l1','l2','l5'],['l2','l4'],['l2','l3'],['l1','l2','l4'],['l1','l3'],['l2','l3'],['l1','l3'],['l1','l2','l3','l5'],['l1','l2','l3']]
		print len(apriori_algorithm(modified_data,min_support))
		#temp=['l1','l2','l3']
		#print (generate_all_subsets(temp,0))
		
	finally:
		data_file.close()
except IOError:
	print('An error occured trying to read the file, file may not be present')



