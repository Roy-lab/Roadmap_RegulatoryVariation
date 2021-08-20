
"""
Defines functions used for running diffusion kernel analysis.
"""
from numpy import *
import scipy as Sci
import scipy.linalg
from scipy.sparse import csr_matrix, dia_matrix, lil_matrix
import time

DEF_TYPE="csr"	# default sparse matrix format

def calculateScores(q, K):
	""" 
	Calculate scores, nothing held aside, of all indices in q.
	Returns a map, score, in format {id : val} 
	For orig, id is an ID, and val is the node's score.
	"""
	score={}
	for v in range(q.shape[0]):
		val = dot(K[v],array(q))
		score[v]=val		
	return score

def crossValidate(val_ids, q, K):
	""" 
	Calculate scores (nothing held aside) and 
	run cross-validation of the list of indices in val_ids.
	Returns 2 maps: (crossval, orig) in format {id : val} 
	For crossval, id is the held-aside ID and val is the node's score when held aside.
	For orig, id is an ID, and val is the node's score in the original.
	"""
	results={}
	orig={}
	for v in range(q.shape[0]):
		if v in val_ids:
			# hold aside v from q
			qprime = array(q)
			qprime[v]=0
			val = dot(K[v],qprime)
			results[v]=val
			oval = val + K[v][v]*q[v]
			orig[v]=oval
			#print v, val, K[v][v], oval
		else:
			#print d
			oval = dot(K[v],array(q))
			orig[v]=oval		
	return (results, orig)

def batchCrossValidate(hitBatches, ordered_names, K):
	""" 
	Given k batches of hits hitBatches= { hitFN : set(hits) }
		- For each batch, remove it from the hit set.
		- calculate scores for all nodes.
	Return a map { hitFN : map}, where each map is { name : val}
	"""
	allHits = set()
	for (fn, hits) in hitBatches.items():
		allHits = set.union(allHits, hits)

	results={}
	for (fn, hits) in hitBatches.items():
		stime = time.time()
		batch = set.difference(allHits, hits)	# hits minus this batch
		q = makeQ(batch, ordered_names)
		resmap = {}
		for v in range(q.shape[0]):
			val = dot(K[v], q)
			resmap[ordered_names[v]]=val
		results[fn]=resmap
		etime = time.time() - stime
		print "# %s processed; %d hits; %.5f seconds" % (fn, len(hits), etime)
	return results

def makeQ(hits, ordered_names):
	"""
	Makes a vertical q vector from the hits.
	For each hit h, the position index(h) in q is 1.
	All other entries are 0.
	"""
	q = []
	for i in range(len(ordered_names)):
		if ordered_names[i] in hits:
			q.append(1)
		else:
			q.append(0)

	return array(q).transpose()
	
def makeQScores(scores, ordered_names):
	"""
	Makes a vertical q vector from an input score map { name : score}.
	For each scored node, the position index(h) in q is score.
	All other entries are 0.
	If a scored node is not present in ordered_names, we skip it.
	"""
	q = []
	for i in range(len(ordered_names)):
		g = ordered_names[i]
		q.append(scores.get(g,0.0))

	return array(q).transpose()
	

def makeKernel(A, lam):
	"""
	Makes a kernel using the symmetric, normalized Laplacian:
	
	L = eye(size(A)) - (D^(-.5))*A*(D^(-.5))
	K = inv(eye(size(L)) + lambda*L)
	
	Input: adjacency matrix (sparse), lambda
	Output: kernel (array))
	"""
	
	# make the ingredients
	D = makeDiagonalDegreeMatrix(A)	
	I = makeSparseEye(A.shape[0])
	isD = invSqrtDiag(D)	# inv sqrt of D

	print "# ... constructed A, D"
	
	# laplacian
	L = I - isD*A*isD	
	print "# ... constructed L. Working on inverting to make K..."

	inside = I + lam*L
	K = linalg.inv(inside.todense())
	return array(K)
	

def makeDiagonalDegreeMatrix(A):
	"""
	Makes a diagonal degree matrix from a sparse adjacency matrix A.
	Returns a s matrix in csr format.
	"""
	# build as lilmatrix
	(rows,cols) = A.nonzero()
	nonzero = zip(rows,cols) # locations of nonzero entries

	D = lil_matrix( (A.shape[0], A.shape[0]))
	# for each row r, place at [r,r] the sum of the columns r,c
	for (r,c) in nonzero:
		D[r,r]+=1
	return D.asformat(DEF_TYPE)

def invSqrtDiag(D):
	""" Makes the square root of D: a diagonal matrix.
	Since it's a diagonal matrix, the inverse sqrt is another
	diagonal matrix. The entries d in D map to (1/sqrt(d)) in the new one.
	"""
	isD = lil_matrix( (D.shape[0], D.shape[0]) )
	(rows,cols) = D.nonzero()
	for (r,c) in zip(rows,cols):
		isD[r,c] = (1.0 / math.sqrt( D[r,c]))
	
	return isD.asformat(DEF_TYPE)
	
def makeSparseMatrix(data):
	"""
	Makes a sparse matrix out of a dictionary.
	Input format:	{ i : {j: w}}
	(where i and j are node indices, and w is a weight)
	
	Output format: matrix in default format)
	"""
	# build as lil_matrix because that's what the wiki told me was fast
	A = lil_matrix( (len(data), len(data)))
	for (i, jmap) in data.items():
		for (j,w) in jmap.items():
			A[i,j]=w
	
	# convert to mtype
	return A.asformat(DEF_TYPE)

def makeSparseEye(size):
	""" Makes a sparse identity matrix of the given size."""
	# make identity matrix
	I = lil_matrix((size,size))
	I.setdiag(ones(size))	
	return I.asformat(DEF_TYPE)

