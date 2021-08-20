#ifndef _EDGE_LIST
#define _EDGE_LIST

#include <map>
#include <vector>
#include <string>

using namespace std;

enum AddType
{
	T_MAX=0,
	T_MIN=1,
	T_SUM=2,
	T_AVG=3
};

class EdgeList
{
public:
	EdgeList(AddType a=T_MAX);
	~EdgeList();
	int addEdge(string tf, string tg, double v);
	int writeToFile(const char*);
private:
	map<string,map<string,double>* > elist;
	map<string,map<string,double>* > ecnt;
	int clear();
	AddType aType;
};

#endif //_EDGE_LIST
