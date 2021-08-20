"""
Given a) a file of gene confidences (including a column indicating which genes 
were used as input to the model that produced the confs)
and b) an external file of genes

produces a list file for use with the AUCPR tool.

Omits input genes from list file.

Oct 2014

Usage:
python make_list_file.py conf_file name_col input_col conf_col gene_list_file name_col2
where conf_file contains at least three columns: 
name_col contains strings
input_col contains floats or ints (nonzero interpreted as input)
conf_col contains floats (confs output by model)
and gene_list_file contains at least one column, name_col2 (strings)
"""
import sys, csv

USAGE="\n".join(("Usage:",
"python make_list_file.py conf_file name_col input_col conf_col gene_list_file name_col2",
"where conf_file contains at least three columns: ",
"name_col contains strings",
"input_col contains floats or ints (nonzero interpreted as input)",
"conf_col contains floats (confs output by model)",
"and gene_list_file contains at least one column, name_col2 (strings)"))

def main(argv):
	if len(argv) != 7:
		print USAGE
		return 2
		
	args=[ a for a in argv ]
	for i in [2,3,4,6]:
		try:
			args[i]=int(args[i])
		except ValueError:
			print USAGE
			return 2			
		
	# confs: { name : {"input":bool, "conf":float}}
	confs=read_map(args[1], args[2], 
	{"input":(lambda x : abs(float(x)) > 0, args[3]), "conf":(lambda x : float(x), args[4])})
	
	# gene set
	genes=read_gene_set(args[5], args[6])
	
	# add missing genes with 0 conf
	for g in genes:
		if g not in confs:
			confs[g]={"input":False, "conf":0.0}
	
	# make the list file contents: (conf, label, name)
	tuples=[]
	for (g, cmap) in confs.items():
		if cmap["input"]:
			continue
		label=0
		if g in genes:
			label=1
		tuples.append((cmap["conf"], label, g))
	
	tuples.sort(key=lambda x : -1*x[0])
	print "#conf\tlabel\tname"
	for t in tuples:
		t1 = "%f\t%d\t%s" % t
		print t1
	
	
def read_map(fn, keycol, col_map):
	"""
	Reads a map from a file. Uses some column (keycol) as key.
	col_map is of the format { "name" : (format_function, column) }
	
	Reads values from other columns into a map, formatted as requested.
 	
	Returns { keycol : { c:row[i] for (c,i) in colmap.items() }}

	"""
	vals={}
	
	with open(fn) as f:
		reader=csv.reader(f, delimiter="\t")	
		
		headerFound=False	
		for row in reader:		
			if row[0][0]=="#":
				continue
				
			name=row[keycol]
			rowmap={}
			for (c, cp) in col_map.items():
				try:
					val = cp[0](row[cp[1]])
					rowmap[c]=val
					
				except ValueError: 
					if headerFound:
						raise Exception("Bad line in file: %s" % row)
						
					headerFound=True
					break	
					
			if len(rowmap) == len(col_map):
				vals[name]=rowmap
#			print name, rowmap.keys()
				
	return vals
	

	

def read_gene_set(fn, col, comment="#", delim="\t"):
	"""
	Reads a gene set from a column of a file.
	Returns set of strings.
	"""
	genes=[]
	with open(fn) as f:
		reader=csv.reader(f, delimiter="\t")	
		for row in reader:
			if row[0][0]=="#":
				continue
			genes.append(row[col])
	return set(genes)
	

if __name__=="__main__":
	sys.exit(main(sys.argv))





