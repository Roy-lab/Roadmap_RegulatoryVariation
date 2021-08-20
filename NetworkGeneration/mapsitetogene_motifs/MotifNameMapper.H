#ifndef _MOTIF_MAPPER
#define _MOTIF_MAPPER

#include <map>
#include <vector>
#include <string>

using namespace std;

class MotifNameMapper
{
public:
	MotifNameMapper();
	~MotifNameMapper();
	int readMap(const char*);
	vector<string>* getTFName(const char*);
private:
	map<string,vector<string>* > mot2tf;
	vector<string>* splitTFs(char* tfnames);
	int clear();
};

#endif //_MOTIF_MAPPER
