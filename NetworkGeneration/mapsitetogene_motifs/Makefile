SRC = Framework.C  EdgeList.C  Gene.C  GeneNameMapper.C  GFFReader.C  MotifNameMapper.C

CC=g++
CFLAGS = -g

mapSiteToGenes: $(SRC)
	$(CC) $(SRC) $(CFLAGS) -o mapSiteToGenes

clean:
	rm -f mapSiteToGenes *~
