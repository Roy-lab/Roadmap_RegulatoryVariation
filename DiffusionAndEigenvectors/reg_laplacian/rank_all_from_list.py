"""
Given an output file with 3 columns (gene, input score, output score), 
and another file with genes to select,
Output and ranks ONLY list genes.
python rank_all_from_list.py scorefile listfile setname
Also indicates whether was input gene.
"""
import sys

if len(sys.argv) != 4:
	print "USAGE: python rank_all_from_list.py scorefile listfile setname"
	sys.exit(2)

# list of genes to focus on
mylist=[]
with open(sys.argv[2]) as f:
	for line in f:
		if line[0]=="#":
			continue
		sp=line.strip().split("\t")
		mylist.append(sp[0])

scores=[]	# [ (gene, score) ]
inputgenes=[]
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
		if inp!=0:
			inputgenes.append(name)
		# skip if not input AND not in list
		if (name not in mylist): # and (inp==0):
			continue
		scores.append( (name, output) )
	
# sort in decreasing order
scores.sort(key=lambda x:-1*x[1])

rank=1
print "name\tscore_%s\trank_%s\tinput_%s" % (sys.argv[3], sys.argv[3],sys.argv[3])
for (gene, score) in scores:
	stat=""
	if gene in inputgenes:
		stat="input"
	print "%s\t%g\t%d\t%s" % (gene,score,rank,stat)
	rank+=1

