#ifndef _NODE_
#define _NODE_
#include <string>
#include <map>
using namespace std;
class Node
{
	public:
	Node();
	~Node();
	int setName(char*);
	int setUnmutatedName(const char*);
	int setTarget(Node* n,double count,int dist, double qvalue);
	int addGenes(char*);
	int setSNP(SNP*);
	int start;
	int end;
	string chrom;
	string name;
	map<string,Node*> targets;
	map<string,double> counts;
	map<string,double> qvals;
	map<string,int> dists;
	map<string,int> geneset;
	map<int,SNP*> snpset;
	private:
};
#endif
