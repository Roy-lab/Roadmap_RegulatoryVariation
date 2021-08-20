"""
Given a file of cross-validation results produced by the reg laplacian kernel code,
produces a list file for use with the AUCPR tool.

Scores for input genes are provided by their leave-out run.
Scores for non-input genes are used from the all-kept run.

Oct 2014

Usage:
python make_list_file_crossval.py 
"""
import sys, csv

USAGE="\n".join(("Usage:",
"python make_list_file_crossval.py file","where file produced by get_kernel_scores.py"))

def main(argv):
	if len(argv) != 2:
		print USAGE
		return 2
		
	# confs: { node : (label, score) }
	confs=read_map(argv[1])
		
	# make the list file contents: (conf, label, name, original_label, original_score)
	tuples = [ (m[0], m[1], name, m[2], m[3]) for (name, m) in confs.items() ]		
	tuples.sort(key=lambda x : -1*x[0])
	print "#conf\tlabel\tname\torig_label\toriginal_score"
	for t in tuples:
		t1 = "%f\t%d\t%s\t%d\t%g" % t
		print t1
	
	
def read_map(fn):
	"""
	Reads map of scores from file produced by get_kernel_scores.py in 
	cross-validation mode. (in addition to columns for name, input score, and 
	output score with none held aside, has one column per hit gene)
	
	##You will need to write a new version of this to handle arbitrary numbers of held-aside genes.## 
	
	Returns a map of confidence scores and labels: { name : {conf, label (1/0)}}
	"""
	vals={}
	
	with open(fn) as f:
		reader=csv.DictReader((row for row in f if not row.startswith('#')), delimiter='\t')	
		
		headers=reader.fieldnames
		# col 0 =gene name
		# col 1 =input score
		# col 2 =none held aside
		# 3+ one gene held aside
		
		# full score column: output_score_LAMBDA
		# leave-out: output_score_LAMBDA_gene (gene may contain "_")
		fullcol=headers[2]				
			
		for row in reader:
			#print row
			# if input gene...
			name=row["name"]
			inscore=float(row["input_score"])
			
			
			# use absolute value as confidence
			fullscore=float(row[fullcol])
			
			# hits could be -1 or 1 -- use absolute value 
			# also use absolute value of held-aside score
			if abs(inscore) > 0.0:
				colname = "%s_%s" % (fullcol, name)
				outscore=float(row[colname])		 ## This won't work if you have multiple genes held aside
				vals[name]=(abs(outscore), 1, inscore, outscore)
				#print name, vals[name], "input"
			else:
				vals[name]=(abs(fullscore), 0, 0, fullscore)
				#print name, vals[name], "not input"
		
				
	return vals
	

	

if __name__=="__main__":
	sys.exit(main(sys.argv))





