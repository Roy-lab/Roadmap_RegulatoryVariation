import sys
import gzip

def readTrans(inname):
	gmap = {}
	f = open(inname,'r')
	for l in f:
		parts = l.strip().split('\t')
		gmap[parts[0]] = parts[1]
	f.close()
	return gmap

def readMotif(inname):
	gmap = {}
	f = open(inname,'r')
	for l in f:
		#M6519_1.02	Tgif2::Tgif1
		parts = l.strip().split('\t')
		gmap[parts[0]] = parts[1].split('::')
	f.close()
	return gmap

def readNet(inname,net,mot2gene,trns2gene):
	f = gzip.open(inname,'r')
	for l in f:
		#M0131_1.02	ENSMUST00000003554	0.00311459
		parts = l.strip().split('\t')
		mot   = parts[0]
		trns  = parts[1]
		v     = float(parts[2])
		if mot in mot2gene and trns in trns2gene:
			g   = trns2gene[trns]
			tfs = mot2gene[mot]
			for tf in tfs:
				if (tf,g) in net:
					if v > net[(tf,g)]:
						net[(tf,g)] = v
				else:
					net[(tf,g)] = v
	f.close()
	return net

def writeNet(outname,net):
	f = open(outname,'w')
	for (tf,tg) in net:
		f.write('%s\t%s\t%f\n' % (tf,tg,net[(tf,tg)]))
	f.close()

if __name__ == '__main__':
	net = {}
	mot2gene  = readMotif(sys.argv[1])
	trns2gene = readTrans(sys.argv[2])
	f = open(sys.argv[3],'r')
	for l in f:
		print l.strip()
		net = readNet(l.strip(),net,mot2gene,trns2gene)
	f.close()
	writeNet(sys.argv[4],net)
