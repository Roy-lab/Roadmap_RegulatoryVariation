#ifndef _FRAMEWORK_
#define _FRAMEWORK_
#include <map>
#include <string>
#include <vector>
#include <fstream>
using namespace std;
class Framework
{
	public:
	Framework();
	~Framework();
	int setMaxDist(int mDist);
	int readPairs(const char*,int,int,int);
	int setRegionNeighborSetSize(int);
	//int createNeighborGraph();
	int generateNeighborGraph(const char* outSuff);
	struct Region
	{
		string chromosome;
		int begin;
		int end;	
	//	int rID;
	};
	struct Pair
	{
		string e;
		string p;
	};
	struct PairIDs
	{
		int e;
		int p;
	};
	private:
	map<string,map<string,double>*> pairSet_True;
	map<string,map<string,double>*> pairSet_Pred;
	map<string,map<string,double>*> pairSet_True_Pro;
	map<string,map<string,double>*> pairSet_Pred_Pro;
	map<string,Region*> regionSet;
	int getDistance(string&,string&);
	map<string,double> geneexp;
	int regionDistNeighborSize;
	//We will use a map with a key per region. For each region we will store a set of regions. The regions that have negative integer
	//keys are before the region and the regions that have positive integer keys are after the region
	map<string,map<string,double>*> neighborSet;
	int maxDist;
	double alpha;
	int maxRegionID;
	map<string,int> allPairs;
	map<int,string> allPairsIDToName;
	map<string,vector<double>*> allPairs_Features;
	map<string,PairIDs*> allPairs_RegionPairs;
	map<int,map<int,int>*> regionPairGraph;
	vector<string> colNames;
};

#endif
