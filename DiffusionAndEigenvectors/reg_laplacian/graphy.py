"""
Contains functions useful for reading/splitting graphs and hit sets.
"""
import time
import sys

def readHits(fnlist, toupper=True):
	"""
	Reads one or more files of hits (provided as list). Returns a map: { name : type }
	where type is one of UP, DOWN, or BOTH.
	if toupper=True, then converts all IDs into uppercase.
	"""
	hits={}
	
	duplicates=set()
	# how many files
	for fn in fnlist:
		fn=fn.strip()
		f = open(fn,"r")
		for line in f:
			if line[0]=="#":
				continue
			sp=line.strip().split("\t")
			name=sp[0]
			if toupper:
				name=name.upper()

			htype=sp[1]
			if name in hits:
				duplicates.add(name)
			hits[name]=htype
		f.close()
	print "# read %d file(s): %s" % (len(fnlist), fnlist)
	print "# found %d duplicate hits" % len(duplicates)
	return hits	

def readHitBatches(fnlist, toupper=True):
	"""
	Reads one or more files of hits (provided as list). Returns a map: { fn : name },
	where fn is the filename.
	"""
	hits={}
	
	for fn in fnlist:
		fn=fn.strip()
		hits[fn]=set()
		f = open(fn,"r")
		for line in f:
			if line[0]=="#":
				continue
			sp=line.strip().split("\t")
			name=sp[0]
			if toupper:
				name=name.upper()

			hits[fn].add(name)
		f.close()
	print "# read %d file(s): %s" % (len(fnlist), fnlist)
	return hits	
	
def readInputScores(fn, cols=[0,1], toupper=True, delim="\t", default=1.0):
	"""
	Reads float scores from tab-delimited file. Returns a map { name : score }.
	Assumes no header line and that all comments start with "#".
	Can provide alternative columns if desired.
	If either only one column provided, then fill with default scores.
	"""
	useDefault=False
	if len(cols)==1:
		useDefault=True
	else:	
		assert(len(cols)==2), "Invalid number of columns provided for readInputScores. Need 1 (for default) or 2 (for scores); provided %d." % len(cols)
	
	scores={}
	with open(fn) as f:
		for line in f:
			line=line.strip()
			if line[0]=="#":
				continue
			sp=line.split(delim)
			name=sp[cols[0]]
			if toupper:	
				name=name.upper()
				
			if useDefault or len(sp)==1:
				val=default
			else:
				val=float(sp[cols[1]])
				
			if name in scores:
				raise Exception("graphy.readInputScores: Duplicate score for item %s." % name)
			scores[name]=val
	return scores
	

def readInteractions(fn, start=0, def_weight=1, toupper=True, read_weight=False, cols=None, delim='\t'):
	""" 
	Reads a file of interactions in the format:
	...\tsource\ttarget\t[weight]
	where "source" is in column start (default: 0)
	may be preceded by an edge ID

	If "cols" provided and has three items, read [col_i, col_j, col_w]

	Returns ( ordered_names, name_map, adj_map )
	where ordered_names is a sorted array of node IDs,
	name_map is a map so we can look up node index by ID,
	and adj_map is a map of the format
	{ i : {j: w}}
	where i is the index into ordered_names, and j is one of i's indices, 
	and w is the weight of the edge.
	***Currently, we treat edges as undirected, and each edge receives weight 1.***
	Ignores self edges!

	"""
	edges=[]
	nodes=set()

	coli = start
	colj = start+1
	colw = start+2

	if (cols != None):
		coli = cols[0]
		colj = cols[1]
		if len(cols) > 2:
			colw = cols[2]
		start=min(cols)
		
	with open(fn) as f:
		for line in f:
			if line[0]=="#":
				continue
			
			# if empty delimiter provided, split on any whitespace
			if len(delim) > 0:
				sp=line.strip().split(delim)
			else:
				sp=line.strip().split()
			if toupper:
				sp=[ s.upper() for s in sp]

			if len(sp) < start+2:
				print >> sys.stderr, "interaction line too short:", sp, fn
				continue

			#for i in sp[start:start+2]:
			#	nodes.add(i)
			i=sp[coli]
			j=sp[colj]
			nodes.add(i)
			nodes.add(j)		
			#edge=(i,j)
			weight=def_weight

			#print i,j

			if len(sp) > colw and read_weight:
				weight=sp[colw]
			edges.append((i,j,weight))

	# create map using names
	adj_map = dict([ (i, {}) for i in nodes])
	for (i,j,w) in edges:
		if i==j:
			continue
		adj_map[i][j]=w
		adj_map[j][i]=w
	#print adj_map
	return adj_map


def readSIF(fn, start=0, def_weight=1, toupper=True, read_weight=False, cols=None, delim="\t"):
	""" 
	Reads a file of interactions in the format: SIF
	...\tsource\tU/D\ttarget
	where "source" is in column start (default: 0)

	Weights are all 1.

	Returns ( ordered_names, name_map, adj_map )
	where ordered_names is a sorted array of node IDs,
	name_map is a map so we can look up node index by ID,
	and adj_map is a map of the format
	{ i : {j: w}}
	where i is the index into ordered_names, and j is one of i's indices, 
	and w is the weight of the edge.
	***Currently, we treat edges as undirected, and each edge receives weight 1.***
	Ignores self edges!

	"""
	edges=[]
	nodes=set()
	f = open(fn,"r")

	coli = start
	colj = start+2

	if (cols != None):
		coli = cols[0]
		colj = cols[1]
		if len(cols) > 2:
			colw = cols[2]
		start=min(cols)

	for line in f:
		if line[0]=="#":
			continue
		sp=line.strip().split(delim)
		if toupper:
			sp=[ s.upper() for s in sp]

		if len(sp) < start+2:
			print >> sys.stderr, "interaction line too short:", sp, fn
			continue

		#for i in sp[start:start+2]:
		#	nodes.add(i)
		i=sp[coli]
		j=sp[colj]
		nodes.add(i)
		nodes.add(j)		
		#edge=(i,j)
		weight=def_weight

		if len(sp) > colw and read_weight:
			weight=sp[colw]
		edges.append((i,j,weight))
	f.close()

	# create map using names
	adj_map = dict([ (i, {}) for i in nodes])
	for (i,j,w) in edges:
		if i==j:
			continue
		adj_map[i][j]=w
		adj_map[j][i]=w
	#print adj_map
	return adj_map

def convertToIndexedMap(adj_map_names):
	"""
	Given an adjacency matrix indexed by names, 
	convert to one indexed by number. 
	Returns:
	ordered_names	ordered list of participants
	name_map		map from name -> ID (reverse of ordered_names)
	adj_map			{ i : { j : w}}, where i and j interact with weight w.			
	"""
	# make name <--> index structures
	ordered_names = sorted(list(adj_map_names.keys()))
	name_map=dict([(ordered_names[i], i) for i in range(len(ordered_names))])

	# create map indexed by number
	# create map using indices
	adj_map = dict([ (i, {}) for i in range(len(ordered_names))])
	for (nameI, imap) in adj_map_names.items():
		for (nameJ, weight) in imap.items():
			(i,j)=(name_map[nameI], name_map[nameJ])
			adj_map[i][j]=weight

	return (ordered_names, name_map, adj_map)

def get_connected_components(adj_map, starts=[], ubiq=[]):
	"""
	Splits an adjacency map into connected components.
	Optional: list of "starts" - nodes we want to cover with our connected
	components. Won't bother searching for components that don't contain any of them.
	Returns a list of adjacency maps.	
	Optional: list of ubiquitous nodes - to ignore during search, but to add back in after split.
	"""
	touched=set()	# items we've touched
	csets = []		# sets of nodes in connected components
	
	# if we specify a hit set, then start from hits
	# otherwise, just go through all nodes in whatever order.
	if len(starts)==0:
		starts=adj_map.keys()

	for n in starts:
		# if we haven't already seen it...
		if n not in touched and n not in ubiq and n in adj_map:
			found = dfs(n, adj_map, ubiq)
			touched = set.union(found, touched)
			csets.append(found)	
	components = []

	for c in csets:
		# get subset of adjacency map
		cmap = dict([ (i, {}) for i in c])
		for i in c:
			for (j, w) in adj_map[i].items():
				if j in cmap:
					cmap[i][j]=w
					cmap[j][i]=w
		components.append(cmap)
	return components

def add_to_map(sub_map, adj_map, to_add):
	""" 
	Given a map, a connected sub-map, and a gene set,
	adds back in connections between the gene set and 
	the connected sub_map.
	Modifies sub_map. Returns the list of added items from to_add.
	"""
	added=[]
	for a in to_add:
		if a not in sub_map:
			sub_map[a]={}
		for (b,w) in adj_map[a].items():
			if b in sub_map:
#				print a, b, w
				sub_map[a][b]=w
				sub_map[b][a]=w
				added.append(a)
	return added


def restrict(adj_map, hits, depth, ubiq=[]):
	""" Returns the adjacency map restricted to nodes within 'depth' hops from any hit."""
	found=set()
	for h in hits:
		if h not in adj_map or h in ubiq:
			continue
		stime =time.time()
		hnodes = dls(h, adj_map, depth, ubiq) 
		etime=time.time()-stime
		print "# \t%s\t%d\t(%.5f s)" % (h, len(hnodes), etime)
		found = set.union(found, hnodes)

	print "%d nodes within distance %d of hits" % (len(found), depth)
	
	new_map={}
	for f in found:
		new_map[f]=dict( [ (j, w) for (j,w) in adj_map[f].items() if j in found])
	return new_map

def dls(start, adj_map, depth, ubiq=[]):
	""" depth limited search"""
	marked=set()
	#print start
	if depth >= 0:
		marked.add(start)
		for j in adj_map[start].keys():
			if j in ubiq:
				continue
			next = dls(j, adj_map, depth-1, ubiq)
			marked = set.union(marked, next)
	#print marked
	return marked

def dfs(start, adj_map, ubiq=[]):
	"""
	dfs to find the connected component that "start" belongs to
	don't include nodes in the ubiquitous list
	"""
	stack = [start]	
	marked = set([start])	

	if start in ubiq:
		return set()

	while len(stack) > 0:		
		v = stack[-1]
		stack=stack[:-1]	
		if v not in adj_map:
			print v	 
		for j in adj_map[v].keys():
			if j in ubiq:
				continue
			if j not in marked :
				marked.add(j)
				stack.append(j)

	return marked	

