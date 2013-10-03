"""
	Usage:
	$python apriori_frequent_itemset.py -f DATASET.csv -s minSupport  -c minConfidence -t template
	
	Eg.
		$ python apriori_frequent_itemset.py -f DATASET.csv -s 0.50 -c 0.70 -t "body has any of G6UP,G82Down,G72UP"

"""


import csv
import copy
import itertools
import re
from optparse 	 import OptionParser

#data_csv_file="GeneAssocaitionData.csv"
gene_prefix="G"
#min_support=0.50
#min_confidence=0.70


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
		

def template1(rules, template):	
	tokens=template.split(" ")
	rule_position=tokens[0]
	quantity=tokens[2]
	items= tokens[4].split(",")
	
	genes=[]
	all_rules=[]
	for item in items:
			genes.append(item)
			
	for rule in rules:	
		body,head=rule
		
		if rule_position=="rule":	
			if quantity=="any":
	
 				combinations=kcombinations(genes,1)
 				for kcombi in combinations:
 					if set(kcombi).issubset(set(body)) or set(kcombi).issubset(set(head)):
 						print rule
 						all_rules.append(rule)
 						
 			elif quantity=="none":
 				is_none=True
 				combinations=kcombinations(genes,1)
				for kcombi in combinations:
 					if not set(kcombi).issubset(set(body)) and set(kcombi).issubset(set(head)):
 						is_none=True
 					else:
 						is_none=False
 						break
 				if is_none:		
 						print rule
 						all_rules.append(rule)
			else:
				k = int(quantity)
 				combinations=kcombinations(genes,k)
 				for kcombi in combinations:
 					if set(kcombi).issubset(set(body)) or set(kcombi).issubset(set(head)):
 						print rule
 						all_rules.append(rule)



		#print body, head
	 	if rule_position=="body":
			if quantity=="any":
 				combinations=kcombinations(genes,1)
 				for kcombi in combinations:
 					if set(kcombi).issubset(set(body)):
 						print rule
 						all_rules.append(rule)
 						
 			elif quantity=="none":
	 			is_none=True
 				combinations=kcombinations(genes,1)
				for kcombi in combinations:
 					if not set(kcombi).issubset(set(body)):
 						is_none=True
 					else:
 						is_none=False
 						break
 				if is_none:		
 						print rule
 						all_rules.append(rule)
			else:
				k = int(quantity)
 				combinations=kcombinations(genes,k)
 				for kcombi in combinations:
 					if set(kcombi).issubset(set(body)):
 						print rule
 						all_rules.append(rule)

 		
 		if rule_position=="head":
			if quantity=="any":
 				combinations=kcombinations(genes,1)
 				for kcombi in combinations:
 					if set(kcombi).issubset(set(head)):
 						print rule
 						all_rules.append(rule)
 						
 			elif quantity=="none":
 				is_none=True
 				combinations=kcombinations(genes,1)
				for kcombi in combinations:
 					if not set(kcombi).issubset(set(head)) :
 						is_none=True 						
 					else:
 						is_none=False
 						break
 				if is_none:		
 						print rule
 						all_rules.append(rule)

			else:
				k = int(quantity)
 				combinations=kcombinations(genes,k)
 				for kcombi in combinations:
 					#print kcombi, head
 					if set(kcombi).issubset(set(head)):
 						print rule
 						all_rules.append(rule)
	
	return all_rules


def template2(rules, template):	
	tokens=template.split(" ")
	sizeOfString=tokens[0]
	startIndex=sizeOfString.find("(")
	endIndex=sizeOfString.find(")")
	
	rule_position=sizeOfString[startIndex+1:endIndex]
	quantity=tokens[2]
	all_rules=[]
	
	for rule in rules:
		body,head=rule
		if rule_position=="rule":
			total_len = len(body)+len(head)
			if total_len >= int(quantity):
				print rule
				all_rules.append(rule)
		
		if rule_position=="body":
			if len(body)>=	int(quantity):
				print rule
				all_rules.append(rule)
		
		if rule_position=="head":
			if len(head)>=	int(quantity):
				print rule
				all_rules.append(rule)
	return all_rules

def template3(rules, template):
		templates=re.split("and|or",template)
		all_rules=[]
		for i, temp in enumerate(templates):
			gen_rules=rules_from_query(rules, temp.strip())
			all_rules.extend(gen_rules)
		
		#and_or=re.search(r'and|or',template,re.M|re.I)
		
		
def rules_from_query(rules, template):
	if "has" in template and not "and" in template and not "or" in template and not "sizeof" in template:
		return template1(rules, template)
	elif "sizeof" in template and not "and" in template and not "or" in template:
		return template2(rules, template)
	else:
		return template3(rules, template)
 
if __name__ == "__main__":

	optparser = OptionParser()
	optparser.add_option('-f', '--inputFile', dest = 'input', help = 'the filename which contains the comma separated values', default="GeneAssocaitionData.csv")
	optparser.add_option('-s', '--minSupport', dest='minS', help = 'minimum support value', default=0.50, type='float')
	optparser.add_option('-c','--minConfidence', dest='minC', help = 'minimum confidence value', default = 0.70, type='float')
	optparser.add_option('-t', '--template', dest = 'template', help = "template to be given in single quotes '....', the genes should be only comma seperated no space allowed", default="body has any of G6UP,G82Down,G72UP")
	
	(options, args) = optparser.parse_args()

	min_support		= options.minS
	min_confidence 	= options.minC				    		
	data_file=open(options.input,"rU")
	line=options.template
	db=read_from_file(data_file)
	rules=(association_rules(db,apriori_algorithm(db,min_support),min_confidence))	
	#line = "body has any of G6UP,G82Down,G72UP";
	all_rules=rules_from_query(rules, line)
	#print_rules(all_rules)


