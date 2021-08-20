"""
Given a set of node scores and a background network (in .sif format),
applies the regularized Laplacian kernel to acquire smoothed
scores for all nodes in the background network. 
Prints tab-delimited scores for all nodes in the background network.

Currently implemented only for UNWEIGHTED graphs, though some functionality
may exist. Don't trust it. :)

The expensive part is calculating the kernel (inverting the matrix).
Rather than doing anything elegant, we just provide the option to save
the kernel object for later (and to read from a pre-created kernel object).

USAGE:
Fresh start: Constructs new kernel from graph.sif and lambda value (> 0), saves the kernel as graph_lambda.npy, and writes smoothed scores to stdout:

python get_kernel_scores.py input_scores.tab graph.sif lambda > output_scores.tab

Read from old kernel: Reads kernel from graph_lambda.npy, writes smoothed scores to stdout. Need to re-supply lambda and graph as a basic compatibility check and in order to match node names properly.

python get_kernel_scores.py input_scores.tab graph.sif lambda graph_lambda.npy > output_scores.tab

Can also provide columns of graph file if desired:

python get_kernel_scores.py input_scores.tab graph.sif --col="0 1" lambda graph_lambda.npy > output_scores.tab

Or standardize capitalization:
--upper

Or run cross-validation: (produces one column per hit)
--crossval

REQUIRES: graphy.py, kernel.py, numpy, scipy

REVISIONS:
Original: chasman@cs.wisc.edu (Sept 2014)
added cross validation (Oct 2014)
Added output filename as last argument (June 2016)
"""
import sys, datetime, os
import numpy as np
import graphy, kernel

USAGE="\n".join(["To construct new kernel: \n\tpython get_kernel_scores.py input_scores.tab graph.sif lambda output_scores.tab", "To read from old kernel: \n\tpython get_kernel_scores.py input_scores.tab graph.sif lambda graph_lambda.npy output_scores.tab", 'Optionally supply graph file columns: --col="0 1"', "Optionally capitalize all identifiers: --upper","Leave-one-out cross-validation on hits: --crossval"])

def main(argv):
	if len(argv) < 5 or len(argv) > 8:
		print >> sys.stderr, USAGE
		return 2
		
	columns=[0,2] # default for sif files
	upperize=False	# check for "--upper" option
	crossval=False	# check for "--crossval" option
	
	# check for columns in argv
	temp=[]
	for a in argv:
		if "--col=" in a:
			sp=a.split("=")
			cs=sp[1].replace('"',"")
			csp=cs.split()
			if len(csp) != 2:
				print >> sys.stderr, "Please supply two columns. Provided %d (command %s)." % (len(csp), a)
				return 2
			columns=[int(x) for x in csp]
		elif "--upper" in a:
			upperize=True
			print "# --upper chosen : Capitalizing all gene names"
		elif "--crossval" in a:
			crossval=True
			print "# --crossval chosen : Running leave-one-out cross-validation on input scores"
		else:
			temp.append(a)		
	argv=temp
			

	# check args to see what we're doing today
	# making a new kernel, or loading an old one
	
	mode="new_k"
	if len(argv) == 5:
		mode="new_k"
	elif len(argv) == 6:
		mode="read_k"

	# check input score file
	scoreFn = argv[1]
	if not os.path.isfile(scoreFn):
		print >> sys.stderr, "Can't find input score file %s." % scoreFn
		return 2
		
	# verify lambda
	lam=1.0
	try:
		lam=float(argv[3])
		
		if lam <= 0:
			raise ValueError
	except ValueError:
		print >> sys.stderr, "Invalid lambda ('%s', converted to '%g'). Please choose lambda > 0.0." % (argv[3], lam)
		return 2
		
		
	# read graph
	graphFn = argv[2]
	if not os.path.isfile(graphFn):
		print >> sys.stderr, "Can't find graph file %s." % graphFn
		return 2
	
	# adjacency map, using names.
	# applies default weight = 1 to all edges!	
	adj_map_names = graphy.readInteractions(graphFn, cols=columns, toupper=upperize, delim="")
	
	# adj map, using numbers
	# ordered_names: maps index to name
	# name_map: name to index 
	# adj_map			{ i : { j : w }}, where i and j interact with weight w.
	# weight functionality not actually implemented into kernel, so don't use it carelessly!
	(ordered_names, name_map, adj_map) = graphy.convertToIndexedMap(adj_map_names)	
		
	# test kernel filename
	ksp = os.path.split(graphFn)	
	testKernFn = "".join(ksp[1].split(".")[:-1]) # without .sif
	
	# if upper chosen, add it to filename
	up=""
	if upperize:
		up="_allcaps"
	testKernFn = os.path.join(ksp[0], "%s_kernel_%g%s.npy" % (testKernFn, lam, up))
	
	# Read input scores (default score if not provided: 1)
	inputScores=None
	try:
		input_scores = graphy.readInputScores(scoreFn, toupper=upperize, default=1.0) # map: hit : score
	except Exception,e:
		print >> sys.stderr, e
		return 1
		
	ingraph = [ h for h in input_scores if h in name_map ]
	print "# Read %d input scores from %s. %d hits in graph." % (len(input_scores), scoreFn, len(ingraph))
	

	# read or create kernel
	if mode=="new_k":
		kernFn = testKernFn					
		print "# Read graph from %s; columns [%s]" % (graphFn, " ".join([str(x) for x in columns]))
		print "# Want to construct kernel with lambda=%g and save it to file %s ... " % (lam, kernFn)
		if os.path.isfile(kernFn):
			print >> sys.stderr, "Kernel object file %s already exists. Delete/move it and re-run to continue." % kernFn
			return 2
								
		# make adjacency matrix	
		A = kernel.makeSparseMatrix(adj_map)
		K = kernel.makeKernel(A, lam)	
		# save
		np.save(kernFn, K)
		print "# Made new kernel (shape %s) and saved to %s" % (str(K.shape), kernFn)
			
	# Read kernel		
	elif mode=="read_k":
		kernFn = argv[4]
		
		# compare graph filename and lambda for compatibility
		if kernFn != testKernFn:
			print >> sys.stderr, "Given your graph filename and lambda, the compatible kernel file should be %s. Instead you provided %s." % (testKernFn, kernFn)
			return 2		
			
		K =	np.load(kernFn)
		print "# Read saved kernel of shape %s from file %s" % (str(K.shape), kernFn)

	
	
	## This is where you would start editing to make 5-fold cross-validation happen
	# Need to construct leaveouts in another way
	# Initialize with an empty set so that we always do the full prediction without leaving out any genes.
	leaveouts=[set()] # Sets of held-aside genes. Default case, none.
					  # Cross-val case, contains sets of individual hits.
	if crossval:
		for h in sorted(input_scores.keys()):
			if h in name_map:
				leaveouts.append(set([h])) 
	#elif n-fold cross validation... (add your code here to construct 5 fold CV)
					  
	scores=[]	      # Same index as leaveouts				  
	for lset in leaveouts:
	
		# make q vector from input scores.
		# first, hold aside any scores.
		new_input = dict( input_scores.items() )
		for i in lset:
			new_input[i]=0
		
		# First, map to indices.		
		q = kernel.makeQScores(new_input, ordered_names)
		
		# intersection between scores and graph
		ingraph = set.intersection(set(name_map.keys()), set(new_input.keys()))
		#print "# %d scored nodes in graph" % (len(ingraph))
	
		# get scores!
		ind_scores = kernel.calculateScores(q, K)
		
		# map scores back to names
		output_scores=dict([ (ordered_names[i], val) for (i,val) in ind_scores.items()])
		scores.append(output_scores)
	
	
	# print scores in node order to output filoe
	outname=argv[-1]
	print "Writing output to %s" % outname
	with open(outname,'w') as f:
	
		# headers: name, input_score, output_score_lamda_(held aside gene(s))
		header=["name", "input_score"]	
		for lset in leaveouts:	
			if len(lset)==0:
				header.append("output_score_%g" % lam)
			else:
				# this might be unwieldy if you have too many to leave out
				header.append("output_score_%g_%s" % (lam, "|".join(lset)))
		print >> f, "\t".join(header)
	
		items=list(ordered_names)
		items.sort()
	
		for name in items:
			slist = [ "%f" % input_scores.get(name, 0.0) ]
			for s in scores:
				slist.append("%f" % s[name] )
			print >> f, "%s\t%s" % (name, "\t".join(slist))


if __name__=="__main__":
	sys.exit(main(sys.argv))

