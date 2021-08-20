import sys

def readMot(inname):
	cmap = {}
	f = open(inname,'r')
	for l in f:
		#chr1	3004298	3004314	0.503	+
		parts = l.strip().split('\t')
		c  = parts[0]
		p1 = int(parts[1])
		p2 = int(parts[2])
		p  = int(p1+p2)/2
		b  = p/1000000
		s  = float(parts[3])
		ss = parts[4]
		t  = cmap.get((c,b),[])
		t.append((p1,p2,s,ss))
		cmap[(c,b)] = t
	f.close()
	return cmap

def readGenes(inname):
	gmap = {}
	f = open(inname,'r')
	for l in f:
		#ENSMUST00000166088	chr10	75032585	75032585
		parts = l.strip().split('\t')
		g  = parts[0]
		c  = parts[1]
		p1 = int(parts[2])
		p2 = int(parts[3])
		b1 = p1/1000000
		b2 = p2/1000000
		gs = gmap.get((c,b1,b2),[])
		gs.append((g,p1,p2))
		gmap[(c,b1,b2)] = gs
	f.close()
	return gmap

def makeNet(cmap,gmap,w,tf):
	net = {}
	for (c,b1,b2) in gmap:
		gs = gmap[(c,b1,b2)]
		for b in range(b1-1,b2+2):
			if (c,b) not in cmap:
				continue
			cs = cmap[(c,b)]
			#print c,b,bb,len(cs),len(gs)
			for (mp1,mp2,v,ss) in cs:
				for (g,p1,p2) in gs:
					#print c,p1,p2,p,g,max(p-w,p1),min(p+w,p2)
					if max(mp1,p1-w) <= min(mp2,p2+w):
						if (tf,g) not in net:
							net[(tf,g)] = v
						else:
							if net[(tf,g)] < v:
								net[(tf,g)] = v
	return net

def writeNet(outname,net):
	f = open(outname,'w')
	for (tf,tg) in net:
		f.write('%s\t%s\t%f\n' % (tf,tg,net[(tf,tg)]))
	f.close()

if __name__ == '__main__':
	cmap = readMot(sys.argv[1])
	gmap = readGenes(sys.argv[2])
	tf   = sys.argv[3]
	w    = int(sys.argv[4])
	net  = makeNet(cmap,gmap,w,tf)
	writeNet(sys.argv[5],net)
