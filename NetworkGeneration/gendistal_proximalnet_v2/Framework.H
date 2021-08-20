#ifndef _FRAMEWORK_
#define _FRAMEWORK_

#include <map>
#include <string>
using namespace std;

class Framework
{
	public:
	Framework();
	~Framework();
	int readInteractions(const char*);
	int readPIQMotifInstances(const char*);
	int generateDistalNetwork(const char*);
	private:
	struct Motif
	{
	string chrom;
	int pos;
	string motifID;
	};

	struct Region
	{
		string chrom;
		int begin;
		int end;
		string name;
	};
	struct Pair
	{
		Region nonpromoter;
		Region promoter;
		string gene;
		string chrom;
	};

	bool getOverlap(Region*, Motif*,string&);
	int getChromInfo(char* tok,string& chr, int& start, int& end);
	map<string,Pair*> interactions;
	map<string,int> enhancerDegree;
	map<string,int> promoterDegree;
	//create a map where the key is the motif name and the entry is a region. 
	map<string,Motif*> sigMotifs;
};

#endif
