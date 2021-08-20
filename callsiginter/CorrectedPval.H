#ifndef _CorrectedPval
#define _CorrectedPval

#include <map>
#include <string>
using namespace std;

struct Pval_T
{
	string name;
	double pval;
	double qval;
	int index;
};

struct customLess
{
	bool operator()(const Pval_T& a, const Pval_T& b) const
	{   
		return a.pval < b.pval;
	}   
};

class CorrectedPval
{
	public:
		CorrectedPval();
		~CorrectedPval();
		int estimateQvalues(map<string,double>& pvals,map<string,double>& qvals);
};

#endif //_CorrectedPval
