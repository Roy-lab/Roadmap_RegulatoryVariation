"""
Given an output file, sorts non-input genes by score.
Outputs only those that did not have input score > 0.

python rank_unknown.py scorefile
"""
import sys

scores=[]	# [ (gene, score) ]
with open(sys.argv[1]) as f:
	for line in f:
		if line[0]=="#":
			continue
		sp=line.strip().split("\t")
		if sp[0]=="name":
			continue
		name=sp[0]
		inp=float(sp[1])
		output=float(sp[2])
		if inp > 0.0 or inp < 0.0:
			continue
		scores.append( (name, output) )
	
# sort in decreasing order
scores.sort(key=lambda x:-1*x[1])

print "name\tscore"
for (gene, score) in scores:
	print "%s\t%g" % (gene,score)

