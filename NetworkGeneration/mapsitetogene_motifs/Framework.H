#ifndef _FRAMEWORK_
#define _FRAMEWORK_
#include <map>
#include <vector>
#include <string>
#include <stdio.h>
using namespace std;
class GFFReader;
class GeneNameMapper;
class MotifNameMapper;
class EdgeList;
class Framework
{
public:
	Framework();
	~Framework();
	int readGFFFile();
	int readSites();
    int mapSitesToGenes();
	int init(int argc, char** argv);
private:
	int printHelp();
	int readRefFasta(const char* aFName,const char*);
	int readChromMap(const char*);
	int readCommonNameMap(const char*);
	int readMotifNameMap(const char*);
	int readSitesSingle(const char*);
	int readSitesDouble(const char*);
	int readSitesBed(const char*);
	struct SiteConf
	{
		int name_index;
		int chr_index;
		int beg_index;
		int end_index;
		int score_index;
		int strand_index;
	} mySiteConf;
	int readBedConfig(const char* aFName);
	struct SiteInfo
	{
		int pos;
		int length;
		double score;
		string chrom;
		char snpClass;
		string name;
		string strand;
	};
	vector<SiteInfo*> siteSet;
	GFFReader gffreader;
	int getChromosomeLength(string&);

	map<string,map<int,int>*> geneSiteMap;
	map<string,int> refchromLength;
	GeneNameMapper gnm;
	MotifNameMapper mnm;
	EdgeList elist;
	map<string,string> srcTargetChromMap;
	char outname[1024];
	char tssname[1024];
	char motname[1024];
	char cnfname[1024];
	char chrmapname[1024];
	char genemapname[1024];
	char motmapname[1024];
	int wsize;
};
#endif
