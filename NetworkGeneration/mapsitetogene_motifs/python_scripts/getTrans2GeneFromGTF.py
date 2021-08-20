import sys
import HTSeq

def getTSS(inname):
	gmap = {}
	gtffile = HTSeq.GFF_Reader(inname)
	for f in gtffile:
		if f.type == "exon" and f.attr["exon_number"] == "1":
			t = f.attr["transcript_id"]
			g = f.attr["gene_name"]
			if t in gmap:
				if gmap[t] != g:
					print 'shouldnt happen!',t,g,gmap[t]
			else:
				gmap[t] = g
	return gmap

def writeGenes(outname,gmap):
	f = open(outname,'w')
	for t in gmap:
		f.write('%s\t%s\n' % (t,gmap[t]))
	f.close()

if __name__ == '__main__':
	gmap = getTSS(sys.argv[1])
	writeGenes(sys.argv[2],gmap)

