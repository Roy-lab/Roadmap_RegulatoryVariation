#ifndef _GFFREADER_
#define _GFFREADER_

#include <map>
#include <string>
using namespace std;

class Gene;

enum TSSType
{
	T_GENE=0,
	T_TRANS=1
};

class GFFReader
{
	public:
	GFFReader(TSSType t=T_GENE);
	~GFFReader();
	int readGFFFile(const char*);
	map<string,Gene*>* getGeneSetForChromosome(string&);
	map<string,map<string,Gene*>*>& getGeneSet();
	int setNeighbors();
	private:
	map<string,map<string,Gene*>*> geneSetPerChrom;
	map<string,Gene*> geneSet;
	map<string,map<int,Gene*>*> geneSetPerChromStrand;
	TSSType tType;
};

#endif
