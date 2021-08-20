#ifndef _GENE_
#define _GENE_
#include <map>
#include <string>
using namespace std;
class Gene
{
	public:
	Gene();
	~Gene();
	int setName(const char*);
	int setGeneCoordinate(int begin, int end);
	int set5utr(int begin, int end);
	int setCDS(int begin,int end);
	int setStrand(char);
	int setChromosome(const char*);
	char getStrand();
	const char* getChromosome();

	int inCDS(int start, int end, int margin);
	int inGene(int start,int end, int margin);
	int in5utr(int start,int end, int margin);
	int get5utrCnts();
	typedef struct
	{
		int begin;
		int end;
	}Coordinate;
	map<string,Coordinate*> cdsSet;
	map<string,Coordinate*> fiveUtrSet;
	int geneStart;
	int geneEnd;
	char strand;
	string chromosome;
	string name;
	Gene* upstreamNeighbor;
	Gene* downstreamNeighbor;
};
#endif
