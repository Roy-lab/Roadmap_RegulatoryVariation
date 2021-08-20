"""
Given a file of cross-validation results produced by the reg laplacian kernel code,
calculates "RMSE" based on the predicted scores vs input scores for held-aside input genes.

Nothing is calculated for genes that were input with 0.

June 2016

Usage:
python calc_rmse.py crossval.tab
"""
import sys, csv, math
from make_list_file_crossval import read_map

USAGE="\n".join(("Usage:",
"python calc_rmse.py crossval_file.tab","where file produced by get_kernel_scores.py"))

def main(argv):
	if len(argv) != 2:
		print USAGE
		return 2
		
	# confs: { node : (score, label, non_abs_label, non_abs_score) }
	confs=read_map(argv[1])
		
	toterr=0.0
	yvals=[]
	for (k, tup) in confs.items():
		# only look at nonzero inputs
		if tup[1]==0:
			continue
		errsq=(tup[2]-tup[3])*(tup[2]-tup[3])
		yvals.append(tup[2])
		#print k, tup[2], tup[3], errsq
		toterr+=errsq
	
	rmse=math.sqrt(toterr)
	yrange=max(yvals)-min(yvals)
	
	# if min==max, then don't actually divide by 0
	if yrange==0:
		yrange=1
	
	
	# normalize by range
	print "%s\t%g\t%g" % (argv[1], rmse/yrange, rmse)
	

if __name__=="__main__":
	sys.exit(main(sys.argv))





