import sys
import HTSeq

def getTSS(inname):
	gmap = {}
	blist = set([])
	gtffile = HTSeq.GFF_Reader(inname)
	for f in gtffile:
		if f.type == "exon" and f.attr["exon_number"] == "1":
			pos = str(f.iv.start_d_as_pos)
			c = 'chr'+pos.split(':')[0]
			p = int(pos.split(':')[1].split('/')[0])
			g = f.attr["transcript_id"]
			if g not in gmap:
				gmap[g] = (c,p,p)
			else:
				(c1,p1,p2) = gmap[g]
				if c != c1:
					print 'shouldnt happen!',c,c1,g,p,p1,p2
					blist.add(g)
					continue
				gmap[g] = (c,min(p1,p),max(p2,p))
	return (gmap,blist)

def writeGenes(outname,gmap,blist):
	f = open(outname,'w')
	for g in gmap:
		if g in blist:
			continue
		(c,p1,p2) = gmap[g]
		f.write('%s\t%s\t%d\t%d\n' % (g,c,p1,p2))
	f.close()

if __name__ == '__main__':
	(gmap,blist) = getTSS(sys.argv[1])
	writeGenes(sys.argv[2],gmap,blist)
