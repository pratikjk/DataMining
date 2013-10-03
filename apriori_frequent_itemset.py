import csv
import copy
import itertools

data_csv_file="GeneAssocaitionData.csv"
gene_prefix="G"
min_support=0.50
min_confidence=0.70


def kcombinations(arr,k):
    """ Returns non empty k combinations of arr"""
    subsets=[]
    for combination in itertools.combinations(arr,k):
    	subsets.append(combination)
    return subsets


def one_item_frequentSet(data_set,min_support):
	""" takes input the data_set and the min_support value to give SET of 1-item_frequent_set"""
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
	""" given two sets and level value, compute all the joins of set1 and set2 of length == level"""
	candidate_set=set()
	for itemset1 in frequent_set1:
		for itemset2 in frequent_set2:
			if len(itemset1.union(itemset2))==level:
				candidate_set.add(itemset1.union(itemset2))
	return candidate_set

def prune(data_set, candidate_set,min_support):
	""" given the entire data_set and the current candidate_set of frequent_items, prune all those items with support less than min_support"""
	pruned_set=set()
	for item_set in candidate_set:
		count=0
		for row in data_set:
			# row is being treated as set here, dunno why
			if item_set.issubset(row):
				count+=1
		if (float(count)/len(data_set)>=min_support):
			pruned_set.add(item_set)
	return pruned_set
		
	

def apriori_algorithm(data_set, min_support):	
	""" given the data_set and the min_support value, find all the frequent items in the dataset """
	L=one_item_frequentSet(data_set,min_support)
	k=1
	temp_candidate_set=set()
	temp_frequent_set=L
	while temp_frequent_set != set([]):
		temp_candidate_set=join(temp_frequent_set,temp_frequent_set,k)
		temp_frequent_set=prune(data_set,temp_candidate_set,min_support)
		L.update(temp_frequent_set)
		k=k+1	
	#print len(L)
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

def prune_for_confidence(data_set, all_subsets, item_set, min_confidence):
	rules=[]
	item_set=set(item_set)
	for subset in all_subsets:
		if len(subset) > 0 and len(subset) < len(item_set):
			subset=set(subset)
			count_head=0
			count_body=0
			body=set(item_set)-set(subset)
			for row in data_set:
				if item_set.issubset(row):
					count_head+=1
				if subset.issubset(row):
					count_body+=1
			if (float(count_head)/float(count_body)) >= min_confidence:
				rule=(tuple(subset),tuple(body))
				rules.append(rule)
	return rules

def association_rules(data_set,all_frequent_item_set,min_confidence):
	# convert the frequent itemsets to list
	temp_set=list(all_frequent_item_set)
	all_frequent_item_list=[]
	for item in temp_set:
    	    all_frequent_item_list.append(list(item))
	
	#print len(all_frequent_item_list)
	#calculate the association rules given the frequent_item_sets
	association_rules=[]
	for item_set in all_frequent_item_list:
			all_subsets=generate_all_subsets(item_set,0)
			rules=prune_for_confidence(data_set,all_subsets,item_set,min_confidence)
			if len(rules)!=0:
				association_rules.extend(rules)
	return association_rules

def read_from_file(data_file):
		reader = csv.reader(data_file)
		modified_data=[]
		for row in reader:
			colnum=1
			modified_row=[]
			for col in row:
				modified_row.append(gene_prefix+str(colnum)+col)
				colnum=colnum+1
			modified_data.append(modified_row)
		return modified_data

def print_rules(rules):
	for rule in rules:
			print "Rule: %s" % (str(rule))
		
def rules_from_query(rules, template):
	
	tokens=template.split(" ")
	rule_position=tokens[0]
	quantity=tokens[2]
	items= tokens[4].split(",")
	
	genes=[]
	for item in items:
			genes.append(item)
	tuple_genes=tuple(genes)
	
	for rule in rules:
		if rule_position=="rule":	
			print 'searching in rule'	
			if quantity=="any":
 				combinations=kcombinations(tuple_genes,1)
 				for kcombi in combinations:
 					if kcombi in rule:
 						print rule
 			elif quantity=="none":
 				combinations=kcombinations(tuple_genes,1)
				for kcombi in combinations:
 					if kcombi not in rule:
 						print rule
			else:
				k = int(quantity)
 				combinations=kcombinations(tuple_genes,k)
 				for kcombi in combinations:
 					if kcombi in rule:
 						print rule


		body,head=rule
		#print body, head
	 	if rule_position=="body":
	 		print 'searching in body'	
			if quantity=="any":
 				combinations=kcombinations(tuple_genes,1)
 				for kcombi in combinations:
 					if kcombi in body:
 						print rule
 			elif quantity=="none":
 				combinations=kcombinations(tuple_genes,1)
				for kcombi in combinations:
 					if kcombi not in body:
 						print rule
			else:
				k = int(quantity)
 				combinations=kcombinations(tuple_genes,k)
 				for kcombi in combinations:
 					if kcombi in body:
 						print rule
 		
 		if rule_position=="head":
	 		print 'searching in head'
			if quantity=="any":
 				combinations=kcombinations(tuple_genes,1)
 				for kcombi in combinations:
 					if kcombi in head:
 						print rule
 			elif quantity=="none":
 				combinations=kcombinations(tuple_genes,1)
				for kcombi in combinations:
 					if kcombi not in head:
 						print rule
			else:
				k = int(quantity)
 				combinations=kcombinations(tuple_genes,k)
 				for kcombi in combinations:
 					if kcombi in head:
 						print rule

# @task: add support for reading filename from commandline input
try:
	data_file=open(data_csv_file,"rU")
	#db=[['l1','l2','l5'],['l2','l4'],['l2','l3'],['l1','l2','l4'],['l1','l3'],['l2','l3'],['l1','l3'],['l1','l2','l3','l5'],['l1','l2','l3']]
	#db=[['M', 'O', 'N', 'K', 'E', 'Y'],['D', 'O', 'N', 'K', 'E', 'Y'],['M', 'A', 'K', 'E'],['M', 'U', 'C', 'K', 'Y'],['C', 'O', 'K', 'I', 'E']]
	db=read_from_file(data_file)
	rules=(association_rules(db,apriori_algorithm(db,min_support),min_confidence))	
	#print_rules(rules)
	line = "head has 1 of G10Down,G72UP";
	rules_from_query(rules, line)
	#kcombinations(('G72UP','G59UP','G96Down'),2)
except IOError:
	print('An error occured trying to read the file, file may not be present')



