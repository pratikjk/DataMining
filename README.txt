README file

Description	: Simple Python implementation of the Apriori Algorithm

Team 		: Pratik Kathalkar
			(pkathalk@buffalo.edu)
			5006-1163
			: Naga Mohan Pokala
			(nagamoha@buffalo.edu)
			5006-1800
Title: Mining Association Rules from Gene Expression Data

Course		: CSE601: DataMining and BioInformatics

Usage:
	$python apriori.py -f DATASET.csv -s minSupport  -c minConfidence -t template
	
	Eg.
		$ python apriori.py -f DATASET.csv -s 0.15 -c 0.6 -t "body has any of G6UP,G82Down,G72UP"


	The program accepts min_support, min_confidence, dataset, and template values through command-line options. If they are not provided, default
	values are also mentioned in the code. The algorithm runs by first generating frequent items-ets.
	These frequent itemsets are then used to generate association rules, by trying out all non-empty subsets of each frequent item sets.
	Those rules which qualify the minimum confidence values are then selected to be the strong rules in the format body==>head.
	
	The template is parsed and then applied on the rules to extract all those rules which satisfy the template conditions.
	Currently only 1st and 2nd template are supported.
	
	All those rules that satisfy the template are printed on the standard output